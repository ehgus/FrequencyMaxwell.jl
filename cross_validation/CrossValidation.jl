"""
# FrequencyMaxwell Cross-Validation Framework

This module provides comprehensive cross-validation capabilities between the Julia 
FrequencyMaxwell implementation and the reference MATLAB ConvergentBornSolver.

The framework ensures numerical accuracy and maintains the 99.9% agreement target 
specified in the project requirements.

## Quick Start

```julia
# Load the framework
include("CrossValidation.jl")
using .CrossValidation

# Set up validation
cv = setup_validation()

# Run quick validation (recommended for development)
summary = run_quick_validation(cv)

# Run comprehensive validation (for final verification)
results = run_comprehensive_validation(cv)

# Validate specific test case
result = validate_single_test(cv, "metalens_forward")

# Check system status
status = get_validation_status()
```

## Framework Structure

- `CrossValidation.jl` - Core validation framework and MATLAB interface
- `ValidationInterface.jl` - High-level user interface functions
- `matlab_scripts/` - MATLAB reference implementations
- `reference_data/` - Saved reference data for offline validation
- `reports/` - Generated validation reports

## Configuration

The framework can be configured through environment variables:

- `MATLAB_PATH` - Path to MATLAB executable
- `MATLAB_SOLVER_PATH` - Path to MATLAB solver directory
- `CV_TOLERANCE` - Default relative tolerance (default: 1e-9)

## Integration with FrequencyMaxwell

The validation framework integrates with the main FrequencyMaxwell package
by comparing results from:

1. Julia implementation using FrequencyMaxwell solvers
2. MATLAB reference implementation using ConvergentBornSolver

Key comparison metrics:
- Field magnitude agreement (>99.9%)
- Phase accuracy (<1e-8 radians)
- Convergence behavior
- Performance benchmarks

## CI/CD Integration

For automated testing in CI/CD pipelines:

```julia
cv = setup_ci_validation()
if cv !== nothing
    results = run_quick_validation(cv, test_subset=:forward)
    @assert results.success_rate >= 0.99 "Validation failed in CI"
end
```
"""
module CrossValidation

# Core functionality
include("src/CrossValidation.jl")
using .CrossValidation

# High-level interface  
include("src/ValidationInterface.jl")
using .ValidationInterface

# Re-export main interface functions
export setup_validation, run_quick_validation, run_comprehensive_validation
export validate_single_test, setup_ci_validation, get_validation_status

# Re-export core types for advanced users
export CrossValidationFramework, ValidationResult, TestCase

# Re-export core validation functions
export run_full_validation, validate_forward_solver, validate_adjoint_solver
export generate_report, load_reference_data, create_test_case

end # module CrossValidation
