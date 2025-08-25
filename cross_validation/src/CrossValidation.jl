"""
Cross-Validation Framework for FrequencyMaxwell Solver

This module provides a comprehensive framework for validating electromagnetic
solvers by comparing results between different implementations.

The framework supports:
- Forward and adjoint solver validation
- Automated testing with configurable parameters
- MATLAB integration via MATLAB.jl
- Reference data generation and comparison
- Detailed validation reporting

GitHub Repository Integration:
This framework automatically downloads and uses the Helmholtz adjoint solver
from: https://github.com/ehgus/Helmholtz-adjoint-solver

Usage:
    using CrossValidation
    
    # Initialize framework (will use GitHub repository)
    cv = CrossValidationFramework()
    
    # Add test cases
    add_test_case!(cv, metalens_forward_test())
    add_test_case!(cv, grating_adjoint_test())
    
    # Run validation
    results = run_validation!(cv)
    
    # Generate report
    generate_report(results, "validation_report.md")

For detailed setup instructions, run: julia setup.jl --help
"""
module CrossValidation

using MATLAB
using JSON3
using Dates
using LinearAlgebra
using Statistics

export CrossValidationFramework, TestCase, ValidationResult
export add_test_case!, run_validation!, generate_report
export metalens_forward_test, grating_forward_test, two_beam_forward_test
export sio2_sphere_forward_test, helical_metalens_forward_test
export metalens_adjoint_test, grating_adjoint_test, double_helix_adjoint_test
export single_helix_adjoint_test

# Load MATLAB integration module
include("MatlabIntegration.jl")
using .MatlabIntegration

# Load enhanced validation metrics module
include("ValidationMetrics.jl")
using .ValidationMetrics

# GitHub repository configuration
const HELMHOLTZ_REPO_URL = "https://github.com/ehgus/Helmholtz-adjoint-solver"
const LOCAL_REPO_PATH = "./Helmholtz-adjoint-solver"

"""
    TestCase

Represents a single validation test case with parameters and expected behavior.

Fields:
- `name::String`: Unique identifier for the test case
- `test_type::String`: Either "forward" or "adjoint"
- `parameters::Dict`: Test-specific parameters
- `tolerances::Dict`: Numerical tolerances for validation
- `description::String`: Human-readable description
"""
struct TestCase
    name::String
    test_type::String  # "forward" or "adjoint"
    parameters::Dict{String, Any}
    tolerances::Dict{String, Float64}
    description::String
end

"""
    ValidationResult

Contains the results of running a validation test case.

Fields:
- `test_case::TestCase`: The original test case
- `success::Bool`: Whether the test passed
- `matlab_results::Dict`: Results from MATLAB solver
- `reference_results::Dict`: Reference/expected results
- `metrics::Dict`: Computed validation metrics
- `execution_time::Float64`: Total execution time in seconds
- `timestamp::DateTime`: When the test was run
- `error_message::String`: Error details if test failed
"""
struct ValidationResult
    test_case::TestCase
    success::Bool
    matlab_results::Dict{String, Any}
    reference_results::Dict{String, Any}
    metrics::Dict{String, Float64}
    execution_time::Float64
    timestamp::DateTime
    error_message::String
end

"""
    CrossValidationFramework

Main framework for managing and executing cross-validation tests.

The framework automatically manages the Helmholtz solver repository and provides
a unified interface for validation testing.

Constructor parameters:
- `matlab_path::String`: Path to MATLAB installation
- `matlab_solver_path::String`: Path to MATLAB solver directory (auto-managed)
- `test_cases::Vector{TestCase}`: List of registered test cases
- `default_tolerances::Dict`: Default numerical tolerances
- `reference_data_path::String`: Path to reference data storage
- `reports_path::String`: Path for validation reports

Example:
    cv = CrossValidationFramework()  # Uses GitHub repository automatically
"""
mutable struct CrossValidationFramework
    matlab_path::String
    matlab_solver_path::String
    test_cases::Vector{TestCase}
    default_tolerances::Dict{String, Float64}
    reference_data_path::String
    reports_path::String
    matlab_session_manager::Union{MatlabSessionManager, Nothing}
    
    function CrossValidationFramework(;
        matlab_path="matlab",
        matlab_solver_path=LOCAL_REPO_PATH,  # Use GitHub repository
        reference_data_path="./reference_data",
        reports_path="./reports"
    )
        
        # Default tolerances for validation
        default_tolerances = Dict(
            "field_relative_error" => 1e-6,
            "field_absolute_error" => 1e-12,
            "energy_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "gradient_relative_error" => 1e-5
        )
        
        new(matlab_path, matlab_solver_path, TestCase[], default_tolerances,
            reference_data_path, reports_path, nothing)
    end
end

"""
    initialize_matlab!(cv::CrossValidationFramework)

Initialize and configure MATLAB session with solver paths using simplified approach.

This function:
1. Creates a simple MATLAB session
2. Adds the solver paths to MATLAB path
3. Verifies that required solver functions are available

Throws an error if MATLAB cannot be initialized or required functions are missing.
"""
function initialize_matlab!(cv::CrossValidationFramework)
    @info "Initializing MATLAB session..."
    
    try
        # Create session manager using simplified approach
        if cv.matlab_session_manager === nothing || !cv.matlab_session_manager.is_active
            cv.matlab_session_manager = initialize_matlab_session()
        end
        
        # Add solver paths to MATLAB path
        solver_path = abspath(cv.matlab_solver_path)
        @info "Adding MATLAB solver path: $solver_path"
        
        if !isdir(solver_path)
            error("Solver directory not found: $solver_path\n" *
                  "Run 'julia setup.jl --clone-repo' to download the repository")
        end
        
        # Add paths recursively
        matlab_session = cv.matlab_session_manager.session
        matlab_cmd = "addpath(genpath('$solver_path'))"
        MATLAB.eval_string(matlab_session, matlab_cmd)
        
        # Also add local matlab_scripts directory
        scripts_path = abspath(joinpath(@__DIR__, "..", "matlab_scripts"))
        if isdir(scripts_path)
            MATLAB.eval_string(matlab_session, "addpath('$scripts_path')")
            @info "Added MATLAB scripts path: $scripts_path"
        end
        
        @info "Solver paths added to MATLAB"
        
        # Verify required functions are available
        required_functions = ["ConvergentBornSolver"]
        available_functions = String[]
        missing_functions = String[]
        
        for func in required_functions
            check_cmd = "exist('$func', 'file')"
            result = MATLAB.eval_string(matlab_session, check_cmd)
            # Convert result to Julia value properly
            if result !== nothing && result > 0
                push!(available_functions, func)
            else
                push!(missing_functions, func)
            end
        end
        
        if !isempty(available_functions)
            @info "Found MATLAB functions: $(join(available_functions, ", "))"
        end
        
        if !isempty(missing_functions)
            @warn "Missing MATLAB functions: $(join(missing_functions, ", "))"
            @warn "Some validation tests may fail. Check solver repository at: $solver_path"
        else
            @info "All required MATLAB functions found"
        end
        
        # Update session manager
        cv.matlab_session_manager.last_used = now()
        
        return true
        
    catch e
        @error "Failed to initialize MATLAB session" exception=e
        @info "For troubleshooting, check that MATLAB.jl is installed: using Pkg; Pkg.add(\"MATLAB\")"
        rethrow(e)
    end
end

"""
    execute_matlab_script(cv::CrossValidationFramework, script_path::String, params::Dict)

Execute a MATLAB script with given parameters through the framework.

Returns the results as a Julia dictionary.
"""
function execute_matlab_script(cv::CrossValidationFramework, script_path::String, params::Dict)
    if cv.matlab_session_manager === nothing || !cv.matlab_session_manager.is_active
        initialize_matlab!(cv)
    end
    
    try
        @info "Executing MATLAB script: $script_path"
        
        # Get active session
        matlab_session = cv.matlab_session_manager.session
        
        # Update last used timestamp
        cv.matlab_session_manager.last_used = now()
        
        # Convert parameters to MATLAB format
        for (key, value) in params
            MATLAB.put_variable(matlab_session, Symbol(key), value)
        end
        
        # Extract script name from path
        script_name = splitext(basename(script_path))[1]
        
        # Build parameter struct in MATLAB
        param_assignments = []
        for (key, value) in params
            if isa(value, String)
                push!(param_assignments, "parameters.$key = '$value';")
            else
                push!(param_assignments, "parameters.$key = $key;")
            end
        end
        
        matlab_code = """
        clear parameters;
        $(join(param_assignments, "\n"))
        
        try
            if exist('$(script_name)', 'file')
                results = $(script_name)('$(params["test_name"])', parameters);
            else
                error('Script $(script_name) not found');
            end
            
            success = true;
            error_msg = '';
        catch ME
            results = struct();
            success = false;
            error_msg = ME.message;
        end
        """
        
        MATLAB.eval_string(matlab_session, matlab_code)
        
        # Get results
        matlab_results = MATLAB.jvalue(MATLAB.get_mvariable(matlab_session, :results))
        success = MATLAB.jvalue(MATLAB.get_mvariable(matlab_session, :success))
        error_msg = MATLAB.jvalue(MATLAB.get_mvariable(matlab_session, :error_msg))
        
        if !success
            error("MATLAB execution failed: $error_msg")
        end
        
        return matlab_results
        
    catch e
        @error "Error executing MATLAB script: $script_path" exception=e
        rethrow(e)
    end
end

"""
    add_test_case!(cv::CrossValidationFramework, test_case::TestCase)

Add a test case to the validation framework.
"""
function add_test_case!(cv::CrossValidationFramework, test_case::TestCase)
    push!(cv.test_cases, test_case)
    @info "Added test case: $(test_case.name) ($(test_case.test_type))"
end

"""
    run_validation!(cv::CrossValidationFramework; 
                   test_filter::Union{String, Nothing}=nothing,
                   save_reference::Bool=false,
                   generate_reports::Bool=true)

Run all registered test cases and return validation results.

Parameters:
- `test_filter`: Run only tests matching this pattern (optional)
- `save_reference`: Save results as new reference data
- `generate_reports`: Automatically generate validation reports

Returns a vector of ValidationResult objects.
"""
function run_validation!(cv::CrossValidationFramework; 
                        test_filter::Union{String, Nothing}=nothing,
                        save_reference::Bool=false,
                        generate_reports::Bool=true)
    
    if isempty(cv.test_cases)
        @warn "No test cases registered. Use add_test_case! to add tests."
        return ValidationResult[]
    end
    
    @info "Starting cross-validation with $(length(cv.test_cases)) test cases"
    
    # Initialize MATLAB if needed
    if cv.matlab_session_manager === nothing || !cv.matlab_session_manager.is_active
        initialize_matlab!(cv)
    end
    
    results = ValidationResult[]
    
    for (i, test_case) in enumerate(cv.test_cases)
        # Apply filter if specified
        if test_filter !== nothing && !occursin(test_filter, test_case.name)
            continue
        end
        
        @info "Running test case $i/$(length(cv.test_cases)): $(test_case.name)"
        
        start_time = time()
        
        try
            # Run the test case
            result = run_single_test(cv, test_case)
            
            execution_time = time() - start_time
            
            # Create validation result
            validation_result = ValidationResult(
                test_case,
                result[:success],
                result[:matlab_results],
                result[:reference_results],
                result[:metrics],
                execution_time,
                now(),
                result[:error_message]
            )
            
            push!(results, validation_result)
            
            status = validation_result.success ? "PASSED" : "FAILED"
            @info "Test $(test_case.name): $status ($(round(execution_time, digits=3))s)"
            
        catch e
            execution_time = time() - start_time
            @error "Test $(test_case.name) failed with exception" exception=e
            
            # Create failed result
            failed_result = ValidationResult(
                test_case,
                false,
                Dict{String, Any}(),
                Dict{String, Any}(),
                Dict{String, Float64}(),
                execution_time,
                now(),
                string(e)
            )
            
            push!(results, failed_result)
        end
    end
    
    @info "Cross-validation completed. $(sum(r.success for r in results))/$(length(results)) tests passed"
    
    # Save reference data if requested
    if save_reference
        save_reference_data(cv, results)
    end
    
    # Generate reports if requested
    if generate_reports
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        report_path = joinpath(cv.reports_path, "validation_report_$timestamp.md")
        generate_report(results, report_path)
    end
    
    return results
end

"""
    run_single_test(cv::CrossValidationFramework, test_case::TestCase)

Execute a single test case and return results.
"""
function run_single_test(cv::CrossValidationFramework, test_case::TestCase)
    try
        # Determine MATLAB script based on test type
        script_name = test_case.test_type == "forward" ? "run_forward_validation" : "run_adjoint_validation"
        script_path = joinpath(cv.matlab_solver_path, "matlab_scripts", "$script_name.m")
        
        # If script doesn't exist in solver directory, use local one
        if !isfile(script_path)
            script_path = joinpath(@__DIR__, "..", "matlab_scripts", "$script_name.m")
        end
        
        # Prepare parameters
        params = copy(test_case.parameters)
        params["test_name"] = test_case.name
        
        # Execute MATLAB script
        matlab_results = execute_matlab_script(cv, script_path, params)
        
        # Load reference data
        reference_results = load_reference_data(cv, test_case.name)
        
        # Compute validation metrics
        metrics = compute_validation_metrics(matlab_results, reference_results, test_case.tolerances)
        
        # Check if validation passed
        success = all(metrics[key] <= tolerance for (key, tolerance) in test_case.tolerances 
                     if haskey(metrics, key))
        
        return Dict(
            :success => success,
            :matlab_results => matlab_results,
            :reference_results => reference_results,
            :metrics => metrics,
            :error_message => success ? "" : "Validation metrics exceeded tolerances"
        )
        
    catch e
        return Dict(
            :success => false,
            :matlab_results => Dict{String, Any}(),
            :reference_results => Dict{String, Any}(),
            :metrics => Dict{String, Float64}(),
            :error_message => string(e)
        )
    end
end

"""
    compute_validation_metrics(matlab_results::Dict, reference_results::Dict, tolerances::Dict)

Compute comprehensive validation metrics by comparing MATLAB results with reference data.
"""
function compute_validation_metrics(matlab_results::Dict, reference_results::Dict, tolerances::Dict)
    try
        # Use enhanced validation metrics (now returns Dict{String, Float64})
        enhanced_metrics = compute_comprehensive_metrics(matlab_results, reference_results, tolerances)
        
        # enhanced_metrics is now already a Dict{String, Float64}
        metrics_dict = Dict{String, Any}(enhanced_metrics)
        
        # Generate quality score for overall assessment
        quality_score = mean([
            get(enhanced_metrics, "numerical_stability", 0.0),
            get(enhanced_metrics, "solution_quality", 0.0), 
            get(enhanced_metrics, "validation_confidence", 0.0)
        ])
        metrics_dict["overall_quality"] = quality_score
        
        return metrics_dict
        
    catch e
        @warn "Error computing enhanced validation metrics, falling back to basic metrics" exception=e
        
        # Fallback to basic metrics computation
        return compute_basic_validation_metrics(matlab_results, reference_results, tolerances)
    end
end

"""
    compute_basic_validation_metrics(matlab_results::Dict, reference_results::Dict, tolerances::Dict)

Compute basic validation metrics as fallback when enhanced metrics fail.
"""
function compute_basic_validation_metrics(matlab_results::Dict, reference_results::Dict, tolerances::Dict)
    metrics = Dict{String, Float64}()
    
    try
        # Basic field comparison
        if haskey(matlab_results, "E_field") && haskey(reference_results, "E_field")
            E_matlab = matlab_results["E_field"]
            E_ref = reference_results["E_field"]
            
            if size(E_matlab) == size(E_ref)
                field_diff = abs.(E_matlab .- E_ref)
                metrics["field_absolute_error"] = maximum(field_diff)
                metrics["field_relative_error"] = maximum(field_diff ./ (abs.(E_ref) .+ 1e-15))
                metrics["field_rms_error"] = sqrt(mean(abs2.(field_diff)))
            end
        end
        
        # Basic energy comparison
        if haskey(matlab_results, "intensity") && haskey(reference_results, "intensity")
            I_matlab = matlab_results["intensity"]
            I_ref = reference_results["intensity"]
            
            if size(I_matlab) == size(I_ref)
                energy_matlab = sum(I_matlab)
                energy_ref = sum(I_ref)
                
                if energy_ref > 0
                    metrics["energy_relative_error"] = abs(energy_matlab - energy_ref) / energy_ref
                end
            end
        end
        
        # Basic convergence metrics
        if haskey(matlab_results, "convergence_history")
            conv_hist = matlab_results["convergence_history"]
            if !isempty(conv_hist)
                metrics["final_residual"] = conv_hist[end]
                metrics["convergence_rate"] = length(conv_hist) > 1 ? 
                    log(conv_hist[1] / conv_hist[end]) / (length(conv_hist) - 1) : 0.0
            end
        end
        
        # Basic efficiency metrics
        for metric_name in ["efficiency", "focusing_efficiency", "scattering_cross_section"]
            if haskey(matlab_results, metric_name) && haskey(reference_results, metric_name)
                val_matlab = matlab_results[metric_name]
                val_ref = reference_results[metric_name]
                
                if val_ref != 0
                    metrics["$(metric_name)_relative_error"] = abs(val_matlab - val_ref) / abs(val_ref)
                end
            end
        end
        
        # Basic gradient validation
        if haskey(matlab_results, "gradient") && haskey(reference_results, "gradient")
            grad_matlab = matlab_results["gradient"]
            grad_ref = reference_results["gradient"]
            
            if size(grad_matlab) == size(grad_ref)
                grad_diff = abs.(grad_matlab .- grad_ref)
                metrics["gradient_absolute_error"] = maximum(grad_diff)
                metrics["gradient_relative_error"] = maximum(grad_diff ./ (abs.(grad_ref) .+ 1e-15))
            end
        end
        
    catch e
        @warn "Error computing basic validation metrics" exception=e
        metrics["computation_error"] = 1.0
    end
    
    return metrics
end

"""
    generate_report(results::Vector{ValidationResult}, output_path::String)

Generate a comprehensive validation report in Markdown format.
"""
function generate_report(results::Vector{ValidationResult}, output_path::String)
    @info "Generating validation report: $output_path"
    
    # Ensure output directory exists
    mkpath(dirname(output_path))
    
    # Generate report content
    report_content = """
# Cross-Validation Report

**Generated:** $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
**GitHub Repository:** $HELMHOLTZ_REPO_URL
**Framework:** FrequencyMaxwell Cross-Validation

## Summary

- **Total Tests:** $(length(results))
- **Passed:** $(sum(r.success for r in results))
- **Failed:** $(sum(!r.success for r in results))
- **Success Rate:** $(round(100 * sum(r.success for r in results) / length(results), digits=1))%

## Test Results

"""
    
    for result in results
        status_emoji = result.success ? "✅" : "❌"
        
        report_content *= """
### $status_emoji $(result.test_case.name)

**Type:** $(result.test_case.test_type)
**Status:** $(result.success ? "PASSED" : "FAILED")
**Execution Time:** $(round(result.execution_time, digits=3))s
**Timestamp:** $(Dates.format(result.timestamp, "yyyy-mm-dd HH:MM:SS"))

**Description:** $(result.test_case.description)

"""
        
        if !isempty(result.metrics)
            report_content *= "**Validation Metrics:**\n"
            for (metric, value) in result.metrics
                tolerance = get(result.test_case.tolerances, metric, NaN)
                status = isnan(tolerance) ? "" : (value <= tolerance ? " ✓" : " ✗")
                report_content *= "- $metric: $(round(value, sigdigits=6))$status\n"
            end
            report_content *= "\n"
        end
        
        if !result.success && !isempty(result.error_message)
            report_content *= "**Error:** $(result.error_message)\n\n"
        end
        
        report_content *= "---\n\n"
    end
    
    # Add configuration information
    report_content *= """
## Configuration

**Repository URL:** $HELMHOLTZ_REPO_URL
**Local Repository Path:** $LOCAL_REPO_PATH
**Test Framework:** Julia CrossValidation.jl
**MATLAB Interface:** MATLAB.jl

## Test Case Details

"""
    
    for result in results
        report_content *= """
### $(result.test_case.name) Parameters

"""
        for (key, value) in result.test_case.parameters
            report_content *= "- **$key:** $value\n"
        end
        
        report_content *= "\n**Tolerances:**\n"
        for (key, value) in result.test_case.tolerances
            report_content *= "- **$key:** $value\n"
        end
        report_content *= "\n"
    end
    
    # Write report to file
    open(output_path, "w") do f
        write(f, report_content)
    end
    
    @info "Report generated successfully: $output_path"
end

"""
    save_reference_data(cv::CrossValidationFramework, results::Vector{ValidationResult})

Save validation results as reference data for future comparisons.
"""
function save_reference_data(cv::CrossValidationFramework, results::Vector{ValidationResult})
    @info "Saving reference data..."
    
    mkpath(cv.reference_data_path)
    
    for result in results
        if result.success
            filename = joinpath(cv.reference_data_path, "$(result.test_case.name)_reference.json")
            
            reference_data = Dict(
                "test_name" => result.test_case.name,
                "test_type" => result.test_case.test_type,
                "timestamp" => string(result.timestamp),
                "matlab_results" => result.matlab_results,
                "parameters" => result.test_case.parameters,
                "repository_url" => HELMHOLTZ_REPO_URL
            )
            
            open(filename, "w") do f
                JSON3.pretty(f, reference_data)
            end
            
            @info "Saved reference data: $(result.test_case.name)"
        end
    end
    
    @info "Reference data saved to $(cv.reference_data_path)"
end

"""
    load_reference_data(cv::CrossValidationFramework, test_case_name::String)

Load reference data for a specific test case.
"""
function load_reference_data(cv::CrossValidationFramework, test_case_name::String)
    filename = joinpath(cv.reference_data_path, "$(test_case_name)_reference.json")
    
    if isfile(filename)
        try
            data = JSON3.read(read(filename, String))
            return data["matlab_results"]
        catch e
            @warn "Failed to load reference data for $test_case_name" exception=e
        end
    end
    
    # Return empty dict if no reference data
    return Dict{String, Any}()
end

# Include test case definitions
include("test_cases.jl")

end # module
