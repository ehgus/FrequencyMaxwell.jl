# FrequencyMaxwell Cross-Validation Framework

A comprehensive framework for validating electromagnetic solvers by comparing results between different implementations. This framework automatically integrates with the **Helmholtz adjoint solver** from the GitHub repository:

**🔗 [https://github.com/ehgus/Helmholtz-adjoint-solver](https://github.com/ehgus/Helmholtz-adjoint-solver)**

## Features

- **Automated GitHub Integration**: Automatically clones and manages the Helmholtz solver repository
- **Forward and Adjoint Solver Validation**: Comprehensive testing of both forward and adjoint electromagnetic solvers
- **MATLAB Integration**: Seamless interface with MATLAB solvers via MATLAB.jl
- **Configurable Test Cases**: Flexible test case definition with customizable parameters and tolerances
- **Automated Reporting**: Generates detailed validation reports in Markdown format
- **Reference Data Management**: Stores and compares against reference solutions
- **Robust Error Handling**: Comprehensive error reporting and troubleshooting guidance

## Quick Start

### Prerequisites

1. **Julia** (>= 1.8) - [Download Julia](https://julialang.org/downloads/)
2. **MATLAB** with valid license - Required for solver execution
3. **Git** - For repository cloning (usually pre-installed on most systems)

### Installation & Setup

1. **Clone this framework** (if not already done):
   ```bash
   git clone <this-repository>
   cd cross_validation
   ```

2. **Setup the framework and download the solver**:
   ```bash
   # Download solver repository and install dependencies
   julia setup.jl --clone-repo --install-deps
   
   # Or run setup interactively
   julia setup.jl
   ```

3. **Run validation tests**:
   ```bash
   # Interactive example with multiple options
   julia example_usage.jl
   
   # Or run directly
   julia -e 'include("example_usage.jl"); run_basic_validation_example()'
   ```

## Repository Structure

```
cross_validation/
├── setup.jl                     # Setup script with GitHub integration
├── example_usage.jl             # Comprehensive usage examples
├── README.md                    # This file
├── src/
│   ├── CrossValidation.jl       # Main framework module
│   ├── ValidationInterface.jl   # Solver interface definitions
│   └── test_cases.jl           # Predefined test case library
├── matlab_scripts/
│   ├── run_forward_validation.m # Forward solver validation scripts
│   └── run_adjoint_validation.m # Adjoint solver validation scripts
├── reference_data/              # Reference solutions (auto-generated)
├── reports/                     # Validation reports (auto-generated)
└── Helmholtz-adjoint-solver/    # GitHub repository (auto-downloaded)
    └── ...                      # Solver implementation
```

## GitHub Repository Integration

The framework automatically manages the [Helmholtz adjoint solver](https://github.com/ehgus/Helmholtz-adjoint-solver) repository:

### Automatic Features

- **Repository Cloning**: Downloads the solver on first use
- **Path Management**: Automatically configures MATLAB paths
- **Update Checking**: Detects when repository updates are available
- **Version Tracking**: Records repository version in validation reports

### Manual Repository Management

```bash
# Clone the repository manually
julia setup.jl --clone-repo

# Check repository status
cd Helmholtz-adjoint-solver
git status
git log --oneline -5

# Update to latest version
git pull origin main

# Switch to specific version/branch
git checkout <branch-name>
git checkout <commit-hash>
```

## Usage Examples

### Basic Usage

```julia
using CrossValidation

# Initialize framework (automatically uses GitHub repository)
cv = CrossValidationFramework()

# Add test cases
add_test_case!(cv, TestCase(
    "metalens_forward_test",
    "forward",
    Dict("wavelength" => 532e-9, "focal_length" => 50e-6),
    Dict("field_relative_error" => 1e-6),
    "Metalens forward solver validation"
))

# Run validation
results = run_validation!(cv)

# Generate report
generate_report(results, "my_validation_report.md")
```

### Advanced Configuration

```julia
# Initialize with custom paths
cv = CrossValidationFramework(
    matlab_solver_path="./Helmholtz-adjoint-solver",  # GitHub repository
    reference_data_path="./my_reference_data",
    reports_path="./my_reports"
)

# Run with filtering and options
results = run_validation!(cv,
    test_filter="metalens",           # Run only tests matching pattern
    save_reference=true,              # Save results as new reference
    generate_reports=true             # Auto-generate reports
)
```

### Available Test Cases

The framework includes pre-defined test cases for common scenarios:

#### Forward Solver Tests
- **Metalens Forward**: Basic metalens focusing validation
- **Grating Forward**: Diffraction grating efficiency validation
- **Two-Beam Forward**: Two-beam interference validation
- **SiO2 Sphere Forward**: Mie scattering validation
- **Helical Metalens Forward**: Vortex beam generation validation

#### Adjoint Solver Tests
- **Metalens Adjoint**: Inverse design optimization validation
- **Grating Adjoint**: Efficiency optimization validation
- **Double Helix Adjoint**: PSF engineering validation
- **Single Helix Adjoint**: Vortex beam optimization validation

## Validation Metrics

The framework computes comprehensive validation metrics:

### Field Validation
- **Field Relative Error**: `max(|E_test - E_ref| / |E_ref|)`
- **Field Absolute Error**: `max(|E_test - E_ref|)`
- **Field RMS Error**: `sqrt(mean(|E_test - E_ref|²))`

### Energy Conservation
- **Energy Relative Error**: `|energy_test - energy_ref| / energy_ref`
- **Intensity Distribution**: Spatial intensity pattern validation

### Convergence Analysis
- **Final Residual**: Final solver residual value
- **Convergence Rate**: Rate of convergence analysis
- **Iteration Count**: Number of iterations to convergence

### Optimization Metrics (Adjoint Tests)
- **Gradient Accuracy**: Adjoint gradient validation
- **Objective Function**: Optimization objective values
- **Design Parameter Updates**: Parameter evolution tracking

## Configuration Options

### Tolerance Settings

```julia
# Custom tolerance configuration
tolerances = Dict{String, Float64}(
    "field_relative_error" => 1e-6,     # Stricter field validation
    "energy_relative_error" => 1e-8,    # Energy conservation
    "convergence_threshold" => 1e-10,   # Solver convergence
    "gradient_relative_error" => 1e-5   # Adjoint gradient accuracy
)
```

### Test Parameters

```julia
# Example test case parameters
parameters = Dict{String, Any}(
    "wavelength" => 532e-9,              # Operating wavelength (m)
    "grid_size" => [64, 64],             # Computational grid size
    "focal_length" => 50e-6,             # Focal length (m)
    "max_iterations" => 100,             # Maximum solver iterations
    "convergence_criterion" => 1e-8      # Solver convergence threshold
)
```

## Repository Integration Details

### Automatic Path Configuration

The framework automatically:
1. Detects the GitHub repository location
2. Adds MATLAB solver paths recursively
3. Verifies required functions are available
4. Provides helpful error messages for missing components

### Repository Requirements

The framework expects the GitHub repository to contain:
- `ConvergentBornSolver.m` or similar solver implementation
- MATLAB function files for forward and adjoint operations
- Standard MATLAB package structure with proper path organization

### Update Management

```bash
# Check for updates
julia setup.jl --check-only

# Update repository
cd Helmholtz-adjoint-solver
git pull origin main

# Verify update
julia setup.jl --test
```

## Troubleshooting

### Common Issues

1. **Repository Not Found**
   ```bash
   # Solution: Download the repository
   julia setup.jl --clone-repo
   ```

2. **MATLAB Not Found**
   ```bash
   # Solution: Install MATLAB and set path
   export MATLAB_PATH=/path/to/matlab
   julia setup.jl --install-deps
   ```

3. **Git Not Available**
   ```bash
   # Solution: Install Git
   # Ubuntu/Debian: sudo apt install git
   # macOS: xcode-select --install
   # Windows: Download from git-scm.com
   ```

4. **Permission Issues**
   ```bash
   # Solution: Fix directory permissions
   chmod +x setup.jl
   chmod -R u+w ./Helmholtz-adjoint-solver
   ```

### Verbose Diagnostics

```bash
# Run setup with detailed output
julia setup.jl --verbose

# Check system requirements
julia setup.jl --check-only --verbose

# Test with verbose output
julia setup.jl --test --verbose
```

### Manual Verification

```julia
# Test framework loading
include("src/CrossValidation.jl")
using .CrossValidation

# Test MATLAB connectivity
using MATLAB
session = MSession()

# Test repository access
isdir("./Helmholtz-adjoint-solver")
readdir("./Helmholtz-adjoint-solver")
```

## Performance Considerations

### Test Execution Time

| Test Type | Grid Size | Typical Time | Memory Usage |
|-----------|-----------|--------------|--------------|
| Forward   | 32×32     | 30-60s       | ~500MB      |
| Forward   | 64×64     | 2-5 min      | ~1-2GB      |
| Adjoint   | 32×32     | 1-3 min      | ~800MB      |
| Adjoint   | 64×64     | 5-15 min     | ~2-4GB      |

### Optimization Tips

1. **Start Small**: Use smaller grid sizes for initial validation
2. **Parallel Execution**: MATLAB can utilize multiple cores automatically  
3. **Memory Management**: Close MATLAB session between large test runs
4. **Incremental Testing**: Use `test_filter` to run specific tests

## Contributing

### Adding New Test Cases

1. Define test case in `src/test_cases.jl`:
   ```julia
   function my_custom_test()
       return TestCase(
           "my_test",
           "forward",  # or "adjoint"
           Dict("parameter" => value),
           Dict("metric" => tolerance),
           "Test description"
       )
   end
   ```

2. Add corresponding MATLAB implementation in `matlab_scripts/`

3. Update exports in `src/CrossValidation.jl`

### Reporting Issues

When reporting issues, please include:
- Julia version (`julia --version`)
- MATLAB version and toolboxes
- Operating system details
- Repository commit hash (`cd Helmholtz-adjoint-solver && git rev-parse HEAD`)
- Complete error messages
- Minimal reproduction example

## Repository Information

- **Primary Repository**: [https://github.com/ehgus/Helmholtz-adjoint-solver](https://github.com/ehgus/Helmholtz-adjoint-solver)
- **Framework Author**: FrequencyMaxwell Cross-Validation Team
- **License**: Follows the license of the Helmholtz solver repository
- **Supported Platforms**: Linux, macOS, Windows (with WSL recommended)

## Related Projects

- [MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) - Julia-MATLAB interface
- [JSON3.jl](https://github.com/quinnj/JSON3.jl) - JSON handling in Julia
- [Julia Language](https://julialang.org/) - High-performance scientific computing

---

For detailed API documentation, run:
```bash
julia --project -e 'using Pkg; Pkg.activate("."); include("src/CrossValidation.jl"); using .CrossValidation; ?CrossValidationFramework'
```

For the latest updates and issues, visit the [Helmholtz adjoint solver repository](https://github.com/ehgus/Helmholtz-adjoint-solver).
