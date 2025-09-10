#!/usr/bin/env julia

"""
Example usage of the simplified cross-validation framework for FrequencyMaxwell.

This demonstrates how to use the simplified MATLAB integration approach
without complex version detection or environment variable management.
"""

# Load the cross-validation framework
include("src/CrossValidation.jl")
using .CrossValidation

# Also load MatlabIntegration for direct access to helper functions
include("src/MatlabIntegration.jl")
using .MatlabIntegration

function main()
    println("="^60)
    println("FrequencyMaxwell Cross-Validation Framework")
    println("Simplified MATLAB Integration Example")
    println("="^60)

    # Check if MATLAB.jl is available
    println("\n1. Checking MATLAB.jl availability...")
    if is_matlab_available()
        println("   ✓ MATLAB.jl is available")
    else
        println("   ✗ MATLAB.jl is not available")
        println("   Install with: using Pkg; Pkg.add(\"MATLAB\")")
        return
    end

    # Create cross-validation framework
    println("\n2. Creating cross-validation framework...")
    cv = CrossValidationFramework()
    println("   ✓ Framework created successfully")

    # Add some test cases
    println("\n3. Adding test cases...")
    test_cases = [
        metalens_forward_test(grid_size = [16, 16], tolerance_level = "relaxed"),
        grating_forward_test(grid_size = [32, 8], tolerance_level = "relaxed"),
        two_beam_forward_test(grid_size = [16, 16], tolerance_level = "relaxed")
    ]

    for test_case in test_cases
        add_test_case!(cv, test_case)
    end

    println("   ✓ Added $(length(cv.test_cases)) test cases")

    # Display configuration
    println("\n4. Framework configuration:")
    println("   MATLAB solver path: $(cv.matlab_solver_path)")
    println("   Reference data path: $(cv.reference_data_path)")
    println("   Reports path: $(cv.reports_path)")

    # Test MATLAB session initialization
    println("\n5. Testing MATLAB session initialization...")
    CrossValidation.initialize_matlab!(cv)
    println("   ✓ MATLAB session initialized successfully")

    # Cleanup
    if cv.matlab_session_manager !== nothing
        cleanup_matlab_session(cv.matlab_session_manager)
        println("   ✓ MATLAB session cleaned up")
    end

    println("\n" * "="^60)
    println("Example completed successfully!")
    println("="^60)
end

# Run the example if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
