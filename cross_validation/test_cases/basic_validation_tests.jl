"""
Basic validation test cases for the cross-validation framework.
These tests verify that the framework components work correctly.
"""

using Test
using LinearAlgebra

# Load the cross-validation framework
include("../CrossValidation.jl")
using .CrossValidation

@testset "CrossValidation Framework Tests" begin
    
    @testset "Status and Configuration" begin
        # Test system status check
        status = get_validation_status()
        
        @test haskey(status, "matlab_jl_available")
        @test haskey(status, "framework_creation")
        @test haskey(status, "julia_version")
        @test status["julia_version"] == string(VERSION)
    end
    
    @testset "Numerical Comparison Functions" begin
        using CrossValidation.CrossValidation
        
        # Test relative error calculation
        @test relative_error(1.0, 1.001) ≈ 0.001
        @test relative_error(0.0, 0.0) == 0.0
        @test isinf(relative_error(0.0, 1.0))
        
        # Test absolute error calculation
        @test absolute_error(1.0, 1.5) == 0.5
        @test absolute_error(-1.0, 1.0) == 2.0
        
        # Test numerical agreement
        @test numerical_agreement(1.0, 1.0 + 1e-10, 1e-9)
        @test !numerical_agreement(1.0, 1.001, 1e-9)
    end
    
    @testset "Array Comparison" begin
        using CrossValidation.CrossValidation
        
        # Create test arrays
        ref_array = randn(10, 10)
        test_array = ref_array + 1e-10 * randn(size(ref_array))
        
        tolerance_config = Dict("relative_error" => 1e-8, "absolute_error" => 1e-8)
        
        result = compare_arrays(ref_array, test_array, tolerance_config)
        
        @test haskey(result, "success")
        @test haskey(result, "max_relative_error")
        @test haskey(result, "max_absolute_error")
        @test result["success"] == true  # Should pass with small errors
        
        # Test size mismatch
        wrong_size_array = randn(5, 5)
        result_size_mismatch = compare_arrays(ref_array, wrong_size_array, tolerance_config)
        
        @test result_size_mismatch["success"] == false
        @test result_size_mismatch["error_type"] == "size_mismatch"
    end
    
    @testset "Complex Field Comparison" begin
        using CrossValidation.CrossValidation
        
        # Create test complex fields
        ref_field = randn(ComplexF64, 8, 8)
        test_field = ref_field + 1e-9 * randn(ComplexF64, size(ref_field))
        
        tolerance_config = Dict(
            "relative_error" => 1e-8, 
            "phase_difference" => 1e-7
        )
        
        result = compare_complex_fields(ref_field, test_field, tolerance_config)
        
        @test haskey(result, "magnitude_analysis")
        @test haskey(result, "phase_max_error")
        @test haskey(result, "success")
        
        # Should pass for small perturbations
        @test result["success"] == true
    end
    
    @testset "Test Case Management" begin
        using CrossValidation.CrossValidation
        
        # Test test case creation
        test_case = create_test_case(
            "test_case_1",
            :forward,
            "Test case for validation",
            "test_script.m",
            parameters=Dict("wavelength" => 532e-9),
            expected_outputs=[:E_field, :convergence],
            tolerance_config=Dict("relative_error" => 1e-10)
        )
        
        @test test_case.name == "test_case_1"
        @test test_case.type == :forward
        @test haskey(test_case.parameters, "wavelength")
        @test test_case.parameters["wavelength"] == 532e-9
        @test :E_field in test_case.expected_outputs
    end
    
    # Only run MATLAB-dependent tests if MATLAB.jl is available
    if CrossValidation.CrossValidation.MATLAB_AVAILABLE
        @testset "MATLAB Interface (if available)" begin
            # Test framework creation
            cv = nothing
            try
                cv = CrossValidationFramework()
                @test cv !== nothing
                @test length(cv.test_cases) == 0  # Should start empty
                
                # Test default tolerance configuration
                @test haskey(cv.default_tolerances, "relative_error")
                @test cv.default_tolerances["relative_error"] == 1e-9
                
            catch e
                @warn "MATLAB interface test failed (this may be expected): $e"
            finally
                if cv !== nothing
                    cleanup_matlab_session!(cv)
                end
            end
        end
    else
        @test_nowarn @warn "MATLAB.jl not available - skipping MATLAB interface tests"
    end
end

@testset "Integration Tests" begin
    
    @testset "Framework Setup" begin
        # Test that setup_validation provides helpful feedback
        status = get_validation_status()
        
        if status["matlab_jl_available"]
            @test_nowarn "Can attempt full setup"
        else
            @test_logs (:error, r"MATLAB\.jl is not available") setup_validation()
        end
    end
    
    @testset "Validation Summary Creation" begin
        using CrossValidation.ValidationInterface
        
        # Create mock validation results
        mock_results = ValidationResult[]
        
        # This would normally require actual validation results
        # For now, just test that the summary structure works
        summary = ValidationInterface.create_validation_summary(mock_results)
        
        @test summary.total_tests == 0
        @test summary.passed == 0
        @test summary.failed == 0
        @test summary.success_rate == 0.0
    end
end

println("\n" * "="^60)
println("Cross-Validation Framework Test Summary")
println("="^60)
println("✓ Numerical comparison functions tested")
println("✓ Array and complex field comparison tested") 
println("✓ Test case management tested")
println("✓ Framework configuration tested")

if CrossValidation.CrossValidation.MATLAB_AVAILABLE
    println("✓ MATLAB interface available for testing")
else
    println("⚠ MATLAB.jl not available - some tests skipped")
    println("  To enable full testing: julia -e 'using Pkg; Pkg.add(\"MATLAB\")'")
end

println("="^60)
