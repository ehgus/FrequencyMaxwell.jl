"""
    ValidationInterface

High-level interface for cross-validation operations.
This module provides simplified functions for common validation tasks.
"""
module ValidationInterface

using ..CrossValidation
using Printf
using Dates
using Statistics

export setup_validation, run_quick_validation, run_comprehensive_validation
export validate_single_test, setup_ci_validation, get_validation_status

"""
    setup_validation(; matlab_path="matlab", force_setup=false)

Set up the cross-validation framework with default configuration.
Returns a configured CrossValidationFramework instance.

# Arguments
- `matlab_path`: Path to MATLAB executable (default: "matlab")
- `force_setup`: Force re-initialization even if already set up

# Examples
```julia
cv = setup_validation()
cv = setup_validation(matlab_path="/usr/local/MATLAB/R2023b/bin/matlab")
```
"""
function setup_validation(; matlab_path = "matlab", force_setup = false)
    @info "Setting up cross-validation framework..."

    # Check if MATLAB is available
    if !is_matlab_available()
        @error """
        MATLAB.jl is not available. To enable cross-validation:
        1. Install MATLAB.jl: julia -e 'using Pkg; Pkg.add("MATLAB")'
        2. Ensure MATLAB is installed and licensed
        3. Restart Julia
        """
        return nothing
    end

    # Create framework instance
    try
        cv = CrossValidationFramework(
            matlab_path = matlab_path,
            matlab_solver_path = abspath(joinpath(@__DIR__, "..", "..", "..",
                "mat-Helmholtz-adjoint-solver")),
            reference_data_path = abspath(joinpath(@__DIR__, "..", "reference_data")),
            reports_path = abspath(joinpath(@__DIR__, "..", "reports"))
        )

        # Initialize MATLAB session
        if !initialize_matlab!(cv)
            @error "Failed to initialize MATLAB session"
            return nothing
        end

        # Load default test cases
        for test_case in create_full_test_suite()
            add_test_case!(cv, test_case)
        end

        @info "Cross-validation framework setup complete!"
        @info "- $(length(cv.test_cases)) test cases registered"
        @info "- MATLAB solver path: $(cv.matlab_solver_path)"
        @info "- Reports will be saved to: $(cv.reports_path)"

        return cv

    catch e
        @error "Failed to setup cross-validation framework" exception=e
        return nothing
    end
end

"""
    run_quick_validation(cv::CrossValidationFramework; test_subset=:forward)

Run a quick validation with a subset of test cases.

# Arguments
- `cv`: Configured CrossValidationFramework
- `test_subset`: Which tests to run (:forward, :adjoint, or :all)

# Returns
- `ValidationSummary`: Summary of results
"""
function run_quick_validation(cv::CrossValidationFramework; test_subset = :forward)
    @info "Running quick validation (subset: $test_subset)"

    # Select test cases
    if test_subset == :forward
        selected_tests = filter(tc -> tc.test_type == "forward", cv.test_cases)[1:min(
            3, end)]
    elseif test_subset == :adjoint
        selected_tests = filter(tc -> tc.test_type == "adjoint", cv.test_cases)[1:min(
            2, end)]
    elseif test_subset == :all
        selected_tests = cv.test_cases[1:min(5, end)]
    else
        error("Invalid test_subset: $test_subset. Use :forward, :adjoint, or :all")
    end

    @info "Selected $(length(selected_tests)) tests for quick validation"

    # Create temporary framework with selected tests
    temp_cv = CrossValidationFramework()
    temp_cv.matlab_session_manager = cv.matlab_session_manager
    for test_case in selected_tests
        add_test_case!(temp_cv, test_case)
    end

    # Run validation
    results = run_validation!(temp_cv, generate_reports = false)

    # Generate summary
    summary = create_validation_summary(results)

    # Quick report
    println("\n" * "="^60)
    println("QUICK VALIDATION SUMMARY")
    println("="^60)
    println("Tests run: $(length(results))")
    println("Passed: $(summary.passed)")
    println("Failed: $(summary.failed)")
    println("Success rate: $(round(summary.success_rate * 100, digits=1))%")
    println("Target (99.9%): $(summary.success_rate >= 0.999 ? "✓ MET" : "✗ NOT MET")")
    println("="^60)

    return summary
end

"""
    run_comprehensive_validation(cv::CrossValidationFramework; generate_report=true)

Run the complete validation suite with all test cases.

# Arguments
- `cv`: Configured CrossValidationFramework  
- `generate_report`: Whether to generate HTML report

# Returns
- `Vector{ValidationResult}`: Complete results
"""
function run_comprehensive_validation(cv::CrossValidationFramework; generate_report = true)
    @info "Running comprehensive validation suite"
    @info "This may take several minutes..."

    start_time = time()

    # Run all validations
    results = run_validation!(cv, generate_reports = false)

    total_time = time() - start_time

    # Generate detailed report
    if generate_report
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
        report_path = joinpath(cv.reports_path, "validation_report_$timestamp.html")

        generate_report(results, report_path)
        @info "Detailed report generated: $report_path"
    end

    # Save reference data for future use
    save_reference_data(cv, results)

    # Print comprehensive summary
    summary = create_validation_summary(results)
    print_comprehensive_summary(summary, total_time)

    return results
end

"""
    validate_single_test(cv::CrossValidationFramework, test_name::String; verbose=true)

Run validation for a single test case with detailed output.
"""
function validate_single_test(cv::CrossValidationFramework, test_name::String; verbose = true)
    # Find test case
    test_case = findfirst(tc -> tc.name == test_name, cv.test_cases)
    if test_case === nothing
        available_tests = [tc.name for tc in cv.test_cases]
        error("Test '$test_name' not found. Available tests: $available_tests")
    end

    @info "Running single test validation: $test_name"

    try
        # Create temporary framework with just this test
        temp_cv = CrossValidationFramework()
        temp_cv.matlab_session_manager = cv.matlab_session_manager
        add_test_case!(temp_cv, cv.test_cases[test_case])

        # Run validation
        results = run_validation!(temp_cv, generate_reports = false)
        result = results[1]  # Get the single result

        if verbose
            print_detailed_test_result(result)
        end

        return result

    catch e
        @error "Test '$test_name' failed" exception=e
        rethrow(e)
    end
end

"""
    setup_ci_validation()

Set up validation for CI/CD environment with minimal dependencies.
"""
function setup_ci_validation()
    @info "Setting up CI/CD validation environment"

    # Check environment
    if !haskey(ENV, "CI") && !haskey(ENV, "GITHUB_ACTIONS")
        @warn "CI environment not detected. This function is designed for automated testing."
    end

    # Attempt setup with timeout
    cv = nothing
    try
        # Set shorter timeout for CI
        cv = setup_validation(force_setup = true)

        if cv === nothing
            @error "CI validation setup failed"
            return nothing
        end

        # Run quick smoke test
        @info "Running CI smoke test..."
        quick_result = run_quick_validation(cv, test_subset = :forward)

        if quick_result.success_rate < 0.5
            @error "CI smoke test failed with success rate: $(quick_result.success_rate)"
            return nothing
        end

        @info "CI validation setup successful"
        return cv

    catch e
        @error "CI setup failed" exception=e
        return nothing
    finally
        # Always cleanup in CI
        if cv !== nothing && cv.matlab_session_manager !== nothing
            cleanup_matlab_session(cv.matlab_session_manager)
        end
    end
end

"""
    get_validation_status()

Get current validation status and system health check.
"""
function get_validation_status()
    status = Dict{String, Any}()

    # Check MATLAB.jl availability
    status["matlab_jl_available"] = is_matlab_available()

    # Check if we can create a framework instance
    cv = CrossValidationFramework()
    status["framework_creation"] = "OK"

    # Try MATLAB session
    if is_matlab_available()
        init_success = initialize_matlab!(cv)
        status["matlab_session"] = init_success ? "OK" : "FAILED"

        if init_success
            cleanup_matlab_session(cv.matlab_session_manager)
        end
    else
        status["matlab_session"] = "UNAVAILABLE"
    end

    # Check paths
    matlab_solver_path = abspath(joinpath(@__DIR__, "..", "..", "..",
        "mat-Helmholtz-adjoint-solver"))
    status["matlab_solver_path"] = isdir(matlab_solver_path) ? "OK" : "MISSING"
    status["solver_path_location"] = matlab_solver_path

    # System info
    status["julia_version"] = string(VERSION)
    status["platform"] = Sys.MACHINE
    status["timestamp"] = string(Dates.now())

    return status
end

# ==================================================================================
# Internal helper functions
# ==================================================================================

struct ValidationSummary
    total_tests::Int
    passed::Int
    failed::Int
    success_rate::Float64
    target_met::Bool
    execution_time::Float64
end

function create_validation_summary(results::Vector{ValidationResult})
    total = length(results)
    passed = count(r -> r.success, results)
    failed = total - passed
    success_rate = total > 0 ? passed / total : 0.0
    target_met = success_rate >= 0.999

    # Average execution time
    avg_time = total > 0 ? mean([r.execution_time for r in results]) : 0.0

    return ValidationSummary(total, passed, failed, success_rate, target_met, avg_time)
end

function print_comprehensive_summary(summary::ValidationSummary, total_time::Float64)
    println("\n" * "="^80)
    println("COMPREHENSIVE VALIDATION SUMMARY")
    println("="^80)
    println(@sprintf("%-25s: %d", "Total tests", summary.total_tests))
    println(@sprintf("%-25s: %d", "Passed", summary.passed))
    println(@sprintf("%-25s: %d", "Failed", summary.failed))
    println(@sprintf("%-25s: %.1f%%", "Success rate", summary.success_rate * 100))
    println(@sprintf("%-25s: %s", "99.9% target",
        summary.target_met ? "✓ MET" : "✗ NOT MET"))
    println(@sprintf("%-25s: %.1f s", "Average test time", summary.execution_time))
    println(@sprintf("%-25s: %.1f s", "Total execution time", total_time))

    if summary.target_met
        println("\n🎉 VALIDATION SUCCESSFUL - All targets met!")
    else
        println("\n⚠️  VALIDATION INCOMPLETE - Review failed tests")
    end
    println("="^80)
end

function print_detailed_test_result(result::ValidationResult)
    println("\n" * "-"^60)
    println("TEST RESULT: $(result.test_case.name)")
    println("-"^60)
    println("Status: $(result.success ? "✓ PASSED" : "✗ FAILED")")
    println("Type: $(result.test_case.test_type)")
    println("Description: $(result.test_case.description)")

    # Timing
    println("\nExecution Time: $(round(result.execution_time, digits=3))s")

    # Metrics
    if !isempty(result.metrics)
        println("\nValidation Metrics:")
        for (metric, value) in result.metrics
            println("  $metric: $(round(value, sigdigits=6))")
        end
    end

    # Error message
    if !result.success && !isempty(result.error_message)
        println("\nError: $(result.error_message)")
    end

    println("-"^60)
end

end # module ValidationInterface
