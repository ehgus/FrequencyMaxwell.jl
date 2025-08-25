"""
Comprehensive Test Suite for Cross-Validation Framework

This test suite validates all components of the cross-validation framework,
including MATLAB integration, validation metrics, error handling, and
end-to-end functionality.
"""

using Test
using Dates
using JSON3

# Add the source directory to load path
push!(LOAD_PATH, "../src")

# Load framework modules
include("../src/CrossValidation.jl")
using .CrossValidation

include("../src/MatlabIntegration.jl")
using .MatlabIntegration

include("../src/ValidationMetrics.jl")
using .ValidationMetrics

include("../src/ErrorAnalysis.jl")
using .ErrorAnalysis

@testset "Cross-Validation Framework Tests" begin
    
    @testset "Test Case Creation" begin
        @test_nowarn metalens_forward_test()
        @test_nowarn grating_forward_test()
        @test_nowarn two_beam_forward_test()
        @test_nowarn sio2_sphere_forward_test()
        @test_nowarn helical_metalens_forward_test()
        
        @test_nowarn metalens_adjoint_test()
        @test_nowarn grating_adjoint_test()
        @test_nowarn double_helix_adjoint_test()
        @test_nowarn single_helix_adjoint_test()
        
        # Test with custom parameters
        custom_test = metalens_forward_test(
            grid_size=[32, 32], 
            wavelength=633e-9, 
            tolerance_level="relaxed"
        )
        @test custom_test.name == "metalens_forward"
        @test custom_test.test_type == "forward"
        @test custom_test.parameters["wavelength"] == 633e-9
        @test custom_test.parameters["grid_size"] == [32, 32]
        
        # Test suite creation
        forward_suite = create_forward_test_suite()
        @test length(forward_suite) == 5
        @test all(tc.test_type == "forward" for tc in forward_suite)
        
        adjoint_suite = create_adjoint_test_suite()
        @test length(adjoint_suite) == 4
        @test all(tc.test_type == "adjoint" for tc in adjoint_suite)
        
        full_suite = create_full_test_suite()
        @test length(full_suite) == 9
        
        quick_suite = create_quick_test_suite()
        @test length(quick_suite) == 3
    end
    
    @testset "CrossValidationFramework Creation" begin
        # Test default constructor
        @test_nowarn cv = CrossValidationFramework()
        
        cv = CrossValidationFramework()
        @test cv.matlab_path == "matlab"
        @test cv.matlab_solver_path == "./Helmholtz-adjoint-solver"
        @test cv.reference_data_path == "./reference_data"
        @test cv.reports_path == "./reports"
        @test isempty(cv.test_cases)
        @test !isempty(cv.default_tolerances)
        
        # Test custom constructor
        custom_cv = CrossValidationFramework(
            matlab_path="/custom/matlab/path",
            reference_data_path="./custom_ref",
            reports_path="./custom_reports"
        )
        @test custom_cv.matlab_path == "/custom/matlab/path"
        @test custom_cv.reference_data_path == "./custom_ref"
        @test custom_cv.reports_path == "./custom_reports"
    end
    
    @testset "Test Case Management" begin
        cv = CrossValidationFramework()
        
        # Test adding test cases
        test_case = metalens_forward_test()
        add_test_case!(cv, test_case)
        @test length(cv.test_cases) == 1
        @test cv.test_cases[1].name == "metalens_forward"
        
        # Add multiple test cases
        add_test_case!(cv, grating_forward_test())
        add_test_case!(cv, metalens_adjoint_test())
        @test length(cv.test_cases) == 3
    end
    
    @testset "MATLAB Integration" begin
        @testset "MATLAB Detection" begin
            # Test MATLAB installation detection
            @test_nowarn detect_matlab_installation()
            
            installation = detect_matlab_installation()
            @test isa(installation, MatlabInstallation)
            @test isa(installation.is_valid, Bool)
            
            # Test validation functions
            if installation.is_valid
                @test !isempty(installation.executable_path)
                @test !isempty(installation.version)
            end
        end
        
        @testset "Session Manager" begin
            # Test session manager creation (without actual MATLAB)
            @test_nowarn MatlabSessionManager()
            
            manager = MatlabSessionManager()
            @test !manager.is_active
            @test isa(manager.session_id, String)
            @test isa(manager.created_at, DateTime)
        end
        
        @testset "Path Detection" begin
            # Test various detection strategies
            @test_nowarn detect_via_environment()
            @test_nowarn detect_via_which_command()
            @test_nowarn detect_via_standard_paths()
            
            # These should return nothing or MatlabInstallation
            env_result = detect_via_environment()
            @test env_result === nothing || isa(env_result, MatlabInstallation)
        end
    end
    
    @testset "Validation Metrics" begin
        @testset "Mock Data Generation" begin
            # Create mock MATLAB results
            matlab_results = Dict{String, Any}(
                "E_field" => complex.(randn(32, 32), randn(32, 32)),
                "intensity" => rand(32, 32),
                "convergence_history" => [1.0, 0.1, 0.01, 0.001, 0.0001],
                "efficiency" => 0.85
            )
            
            # Create mock reference results with small differences
            reference_results = Dict{String, Any}(
                "E_field" => matlab_results["E_field"] + 1e-8 * complex.(randn(32, 32), randn(32, 32)),
                "intensity" => matlab_results["intensity"] + 1e-8 * randn(32, 32),
                "efficiency" => 0.851
            )
            
            tolerances = Dict{String, Float64}(
                "field_relative_error" => 1e-6,
                "energy_relative_error" => 1e-6
            )
            
            @test_nowarn metrics = compute_comprehensive_metrics(
                matlab_results, reference_results, tolerances
            )
            
            metrics = compute_comprehensive_metrics(matlab_results, reference_results, tolerances)
            @test isa(metrics, ValidationMetrics)
            
            # Test individual metric computation functions
            @test_nowarn compute_field_metrics(matlab_results, reference_results)
            @test_nowarn compute_energy_metrics(matlab_results, reference_results)
            @test_nowarn compute_convergence_metrics(matlab_results)
            @test_nowarn compute_statistical_metrics(matlab_results, reference_results)
        end
        
        @testset "Convergence Analysis" begin
            # Test convergence analysis
            conv_history = [1.0, 0.5, 0.1, 0.05, 0.01, 0.001]
            analysis = analyze_convergence(conv_history)
            
            @test haskey(analysis, "initial_residual")
            @test haskey(analysis, "final_residual")
            @test haskey(analysis, "reduction_factor")
            @test analysis["initial_residual"] == 1.0
            @test analysis["final_residual"] == 0.001
            @test analysis["reduction_factor"] ≈ 1000.0
        end
        
        @testset "Energy Conservation" begin
            # Test energy conservation checking
            results_with_energy = Dict{String, Any}(
                "incident_power" => 1.0,
                "scattered_power" => 0.7,
                "absorbed_power" => 0.3
            )
            
            diagnostics = check_energy_conservation(results_with_energy)
            @test haskey(diagnostics, "conservation_error")
            @test haskey(diagnostics, "conserved")
            @test diagnostics["conserved"] == true  # Perfect conservation
        end
    end
    
    @testset "Error Analysis" begin
        @testset "Error Classification" begin
            # Test error classification
            error_type, category, severity = classify_error("MATLAB not found", Dict{String, Any}())
            @test error_type == "MATLAB_NOT_FOUND"
            @test category == "MATLAB_INSTALLATION"
            @test severity == "CRITICAL"
            
            error_type, category, severity = classify_error("validation tolerance exceeded", Dict{String, Any}())
            @test error_type == "VALIDATION_TOLERANCE"
            @test category == "NUMERICAL_ACCURACY"
            @test severity == "MEDIUM"
        end
        
        @testset "Error Report Creation" begin
            # Test error report creation
            test_error = ValidationError("Test validation failed")
            context = Dict{String, Any}("test_name" => "test_case")
            
            @test_nowarn ErrorReport(test_error, context)
            
            error_report = ErrorReport(test_error, context)
            @test isa(error_report, ErrorReport)
            @test !isempty(error_report.error_message)
            @test !isempty(error_report.error_type)
            @test isa(error_report.timestamp, DateTime)
        end
        
        @testset "System Diagnostics" begin
            # Test system information collection
            @test_nowarn collect_system_information()
            
            sys_info = collect_system_information()
            @test haskey(sys_info, "os")
            @test haskey(sys_info, "architecture")
            @test haskey(sys_info, "julia_version")
            @test haskey(sys_info, "cpu_cores")
        end
    end
    
    @testset "Validation Result Structure" begin
        # Test ValidationResult creation
        test_case = metalens_forward_test()
        matlab_results = Dict{String, Any}("E_field" => rand(16, 16))
        reference_results = Dict{String, Any}("E_field" => rand(16, 16))
        metrics = Dict{String, Float64}("field_relative_error" => 1e-6)
        
        result = ValidationResult(
            test_case,
            true,
            matlab_results,
            reference_results,
            metrics,
            1.5,
            now(),
            ""
        )
        
        @test result.success == true
        @test result.execution_time == 1.5
        @test result.test_case.name == "metalens_forward"
        @test !isempty(result.matlab_results)
    end
    
    @testset "File I/O Operations" begin
        # Create temporary directories for testing
        temp_ref_dir = mktempdir()
        temp_reports_dir = mktempdir()
        
        try
            cv = CrossValidationFramework(
                reference_data_path=temp_ref_dir,
                reports_path=temp_reports_dir
            )
            
            # Test reference data saving/loading
            test_case = metalens_forward_test()
            matlab_results = Dict{String, Any}(
                "E_field" => rand(8, 8),
                "intensity" => rand(8, 8),
                "efficiency" => 0.8
            )
            
            result = ValidationResult(
                test_case, true, matlab_results, Dict{String, Any}(), 
                Dict{String, Float64}(), 1.0, now(), ""
            )
            
            # Test saving reference data
            @test_nowarn save_reference_data(cv, [result])
            
            # Verify file was created
            ref_file = joinpath(temp_ref_dir, "metalens_forward_reference.json")
            @test isfile(ref_file)
            
            # Test loading reference data
            @test_nowarn loaded_ref = load_reference_data(cv, "metalens_forward")
            loaded_ref = load_reference_data(cv, "metalens_forward")
            @test haskey(loaded_ref, "E_field")
            
            # Test report generation
            @test_nowarn generate_report([result], joinpath(temp_reports_dir, "test_report.md"))
            
            report_file = joinpath(temp_reports_dir, "test_report.md")
            @test isfile(report_file)
            
            # Verify report content
            report_content = read(report_file, String)
            @test occursin("Cross-Validation Report", report_content)
            @test occursin("metalens_forward", report_content)
            
        finally
            # Cleanup
            rm(temp_ref_dir, recursive=true, force=true)
            rm(temp_reports_dir, recursive=true, force=true)
        end
    end
    
    @testset "Framework Configuration" begin
        # Test default tolerances
        cv = CrossValidationFramework()
        @test haskey(cv.default_tolerances, "field_relative_error")
        @test haskey(cv.default_tolerances, "energy_relative_error")
        @test cv.default_tolerances["field_relative_error"] == 1e-6
        
        # Test tolerance validation
        test_case = metalens_forward_test(tolerance_level="strict")
        @test test_case.tolerances["field_relative_error"] == 1e-8
        
        test_case_relaxed = metalens_forward_test(tolerance_level="relaxed")
        @test test_case_relaxed.tolerances["field_relative_error"] == 1e-3
    end
    
    @testset "Mock Validation Run" begin
        # Test a complete mock validation run (without MATLAB)
        cv = CrossValidationFramework()
        
        # Add a quick test case
        add_test_case!(cv, metalens_forward_test(grid_size=[8, 8], tolerance_level="relaxed"))
        
        @test length(cv.test_cases) == 1
        
        # This would normally call run_validation!, but we can't without MATLAB
        # Instead, test the internal components that would be used
        test_case = cv.test_cases[1]
        @test test_case.name == "metalens_forward"
        @test test_case.test_type == "forward"
        @test test_case.parameters["grid_size"] == [8, 8]
    end
    
    @testset "Error Handling" begin
        # Test error handling in various scenarios
        cv = CrossValidationFramework()
        
        # Test empty test case list
        results = run_validation!(cv)
        @test isempty(results)
        
        # Test invalid reference data loading
        invalid_ref = load_reference_data(cv, "nonexistent_test")
        @test isempty(invalid_ref)
        
        # Test metric computation with invalid data
        empty_matlab = Dict{String, Any}()
        empty_ref = Dict{String, Any}()
        tolerances = Dict{String, Float64}()
        
        metrics = compute_basic_validation_metrics(empty_matlab, empty_ref, tolerances)
        @test isa(metrics, Dict)
    end
    
    @testset "Integration with Test Cases" begin
        # Test that all predefined test cases are valid
        test_constructors = [
            metalens_forward_test,
            grating_forward_test,
            two_beam_forward_test,
            sio2_sphere_forward_test,
            helical_metalens_forward_test,
            metalens_adjoint_test,
            grating_adjoint_test,
            double_helix_adjoint_test,
            single_helix_adjoint_test
        ]
        
        for constructor in test_constructors
            @test_nowarn test_case = constructor()
            test_case = constructor()
            @test !isempty(test_case.name)
            @test test_case.test_type in ["forward", "adjoint"]
            @test !isempty(test_case.parameters)
            @test !isempty(test_case.tolerances)
            @test !isempty(test_case.description)
        end
    end
end

# Additional integration tests that require manual verification
@testset "Integration Tests (Manual)" begin
    @testset "Repository Verification" begin
        # These tests check repository integration but don't require MATLAB
        repo_path = "./Helmholtz-adjoint-solver"
        
        if isdir(repo_path)
            @test isdir(joinpath(repo_path, ".git"))
            
            # Check for key files
            key_files = ["README.md", "LICENSE"]
            for file in key_files
                found = false
                for (root, dirs, files) in walkdir(repo_path)
                    if file in files
                        found = true
                        break
                    end
                end
                if !found
                    @info "Key file $file not found in repository"
                end
            end
        else
            @info "Repository not found at $repo_path - run 'julia setup.jl --clone-repo'"
        end
    end
    
    @testset "System Requirements Check" begin
        # Check Julia version
        @test VERSION >= v"1.8"
        
        # Check available memory (basic check)
        if Sys.total_memory() < 2^32  # 4 GB
            @info "System has less than 4 GB RAM - some tests may fail"
        end
        
        # Check disk space (basic check)
        try
            if Sys.isunix()
                df_output = read(`df -B1 .`, String)
                @test !isempty(df_output)
            end
        catch
            @info "Could not check disk space"
        end
    end
end

println("Cross-validation framework test suite completed.")
println("Run with MATLAB available to test full functionality.")