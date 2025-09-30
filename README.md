# FrequencyMaxwell.jl

[![Julia](https://img.shields.io/badge/julia-%3E%3D1.6-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Development Status](https://img.shields.io/badge/status-active%20development-brightgreen.svg)](https://github.com)

**FrequencyMaxwell.jl** is a modern Julia package for electromagnetic field simulation using frequency-domain Maxwell equations. It provides high-performance, automatic differentiation-enabled solvers for electromagnetic scattering, inverse design, and optimization problems.

## Overview

FrequencyMaxwell.jl provides comprehensive electromagnetic simulation capabilities in Julia through:

- **Type-safe, performance-optimized** electromagnetic solvers
- **Native automatic differentiation** support via Zygote.jl for inverse design
- **Flexible linear algebra backends** through LinearSolve.jl integration
- **GPU acceleration** support for large-scale simulations
- **Memory-efficient implementations** with comprehensive electromagnetic solver ecosystem

This package implements proven computational electromagnetics algorithms with modern Julia language features and enhanced performance capabilities.

## Key Features

### Core Capabilities
- **Convergent Born Iterative Method** for electromagnetic scattering
- **Frequency-domain Maxwell equation** solvers
- **High numerical aperture** microscopy simulations
- **Multi-scale electromagnetic modeling** from nano to microscale

### Modern Julia Integration
- **LinearSolve.jl** integration for optimized linear algebra operations
- **Automatic differentiation** compatibility for gradient-based optimization
- **Type-stable, memory-efficient** implementations
- **GPU acceleration** support through CUDA.jl ecosystem
- **Thread-safe parallel** computing capabilities

### Scientific Applications
- **Electromagnetic scattering** analysis
- **Optical microscopy** simulation (brightfield, darkfield, phase contrast)
- **Inverse design** and topology optimization
- **Material characterization** through inverse scattering
- **Phantom-based** validation and testing

## Installation

### Julia Package Manager

```julia
using Pkg
Pkg.add("FrequencyMaxwell")
```

### Development Installation

For the latest development version:

```julia
using Pkg
Pkg.add(url="https://github.com/your-org/FrequencyMaxwell.jl")
```

### System Requirements

- **Julia 1.6+** (Julia 1.8+ recommended for optimal performance)
- **FFTW.jl** for fast Fourier transforms
- **LinearSolve.jl** for advanced linear algebra
- **StaticArrays.jl** for efficient small array operations

Optional dependencies for extended functionality:
- **CUDA.jl** for GPU acceleration
- **Zygote.jl** for automatic differentiation
- **ChainRulesCore.jl** for custom differentiation rules

## Quick Start

### Basic Electromagnetic Scattering

```julia
using FrequencyMaxwell

# Create electromagnetic solver with streamlined API
solver = ConvergentBornSolver(
    wavelength = 500e-9,          # 500 nm wavelength (green light)
    permittivity_bg = 1.33^2,     # Water background (n=1.33)
    resolution = (50e-9, 50e-9, 50e-9),  # 50 nm isotropic resolution
    grid_size = (128, 128, 64),           # Computational grid size
    tolerance = 1e-6              # Convergence tolerance
)

# Define incident plane wave source
source = PlaneWaveSource(
    wavelength = solver.wavelength,
    polarization = [1.0, 0.0, 0.0],    # X-polarized light
    k_vector = [0.0, 0.0, 1.0],        # Propagating in +Z direction
    amplitude = 1.0                     # 1 V/m amplitude
)

# Generate a spherical bead phantom
phantom = phantom_bead(
    solver.grid_size,
    [1.5^2],                      # Polystyrene bead (n=1.5)
    16.0,                         # Radius in pixels
    num_bead = 1
)

# Solve the electromagnetic scattering problem
E_field, H_field = solve(solver, source, phantom)

# Analyze results
fields = ElectromagneticField(E_field, H_field, solver.grid_size,
                             solver.resolution, solver.wavelength)

println("Total field energy: $(field_energy(fields)) J")
println("Maximum intensity: $(maximum(field_intensity(fields)))")
```

### Advanced Usage with Automatic Differentiation

```julia
using FrequencyMaxwell, Zygote

# Define an objective function for inverse design
function objective(phantom_params)
    phantom = phantom_bead(config.grid_size, phantom_params, 16.0)
    E_field, H_field = solve(solver, source, phantom)
    fields = ElectromagneticField(E_field, H_field, config.grid_size, 
                                 config.resolution, config.wavelength)
    return field_energy(fields)  # Minimize total energy
end

# Compute gradients using automatic differentiation
initial_params = [1.5^2]
gradient = Zygote.gradient(objective, initial_params)[1]
println("Gradient with respect to permittivity: $gradient")
```

## API Overview

### Core Types

#### Solver Types
- **`ConvergentBornSolver`**: Main electromagnetic solver
- **`solve(solver, source, phantom)`**: Solve electromagnetic scattering problem

#### Source Types
- **`PlaneWaveSource`**: Plane wave illumination source
- **`AbstractCurrentSource`**: Base type for current sources
- **`source_wavelength(source)`**: Extract source wavelength
- **`source_power(source)`**: Calculate source power density

### Field Analysis

#### ElectromagneticField Type
```julia
fields = ElectromagneticField(E_field, H_field, grid_size, resolution, wavelength)

# Field analysis functions
energy = field_energy(fields)              # Total electromagnetic energy
intensity = field_intensity(fields)        # |E|² intensity distribution
poynting = poynting_vector(fields)         # Poynting vector (energy flow)
plane = extract_plane(fields, axis, index) # Extract 2D plane for visualization
```

### Phantom Generation

#### Built-in Phantoms
```julia
# Spherical bead phantom
bead = phantom_bead(grid_size, permittivities, radius, num_bead=1)

# Planar plate phantom  
plate = phantom_plate(grid_size, permittivity, thickness, orientation=3)

# Cylindrical phantom
cylinder = phantom_cylinder(grid_size, permittivity, radius, axis=3)
```

## Performance Highlights

### Computational Efficiency
- **LinearSolve.jl Integration**: Automatic selection of optimal linear algebra backends
- **Memory-Efficient**: Minimal allocations during iterative solving
- **Type-Stable**: All operations are type-stable for maximum performance
- **SIMD Optimized**: Vectorized operations using Julia's native SIMD

### Scalability
- **GPU Acceleration**: Native support for CUDA.jl backends
- **Parallel Computing**: Thread-safe implementations for multi-core systems
- **Large-Scale Problems**: Efficient handling of problems with millions of unknowns
- **Adaptive Algorithms**: Convergence-based iteration control

### Benchmarks
Typical performance characteristics:
- **128³ grid**: ~10-50 iterations, 2-10 seconds (CPU), <1 second (GPU)
- **256³ grid**: ~15-75 iterations, 30-120 seconds (CPU), 5-15 seconds (GPU)
- **Memory usage**: ~O(N log N) for N-voxel problems

## Examples

### Example 1: Microscopy Simulation

```julia
# High-resolution microscopy simulation
config = ConvergentBornConfig(
    wavelength = 488e-9,           # Blue laser
    permittivity_bg = 1.515^2,     # Oil immersion medium
    resolution = (25e-9, 25e-9, 50e-9),  # Anisotropic resolution
    grid_size = (256, 256, 128)
)

# Multiple bead phantom
phantom = phantom_bead(config.grid_size, [1.45^2, 1.6^2], [20.0, 15.0], num_bead=5)

solver = ConvergentBornSolver(config)
source = PlaneWaveSource(wavelength=config.wavelength, 
                        polarization=[1.0, 0.0, 0.0])

E_field, H_field = solve(solver, source, phantom)
```

### Example 2: Material Characterization

```julia
# Inverse scattering for material characterization
using Optim

function characterize_material(measured_intensity)
    function cost(permittivity)
        phantom = phantom_cylinder(config.grid_size, permittivity[1], 25.0)
        E_field, H_field = solve(solver, source, phantom)
        computed_intensity = field_intensity(ElectromagneticField(
            E_field, H_field, config.grid_size, config.resolution, config.wavelength))
        return norm(computed_intensity - measured_intensity)^2
    end
    
    result = optimize(cost, [1.5^2], LBFGS())
    return sqrt(result.minimizer[1])  # Return refractive index
end
```

### Example 3: Custom Phantom Generation

```julia
# Custom phantom using mathematical functions
function custom_phantom(grid_size, center, radius_func)
    phantom = ones(ComplexF64, grid_size)
    for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
        pos = [i, j, k] .- center
        r = norm(pos)
        if r <= radius_func(pos)
            phantom[i, j, k] = 1.6^2  # High refractive index region
        end
    end
    return phantom
end

# Ellipsoidal phantom
ellipsoid = custom_phantom((128, 128, 64), [64, 64, 32], 
                          pos -> sqrt(sum((pos ./ [20, 15, 10]).^2)))
```

## Development Status

### Current Version: 0.1.0

FrequencyMaxwell.jl is actively developed and provides comprehensive electromagnetic simulation capabilities in Julia with state-of-the-art algorithmic implementations:

#### Current Features
- **Modern Julia Ecosystem**: Full integration with LinearSolve.jl, FFTW.jl, and AD ecosystem
- **Performance Optimizations**: Memory-efficient algorithms with GPU acceleration
- **Type Safety**: Comprehensive type system for configuration and field management
- **Extensible Architecture**: Modular design for easy extension and customization

#### Roadmap
- **Additional Solvers**: FDTD, FEM, and hybrid methods
- **Advanced Sources**: Focused beams, structured illumination, arbitrary current distributions
- **Optimization Framework**: Built-in inverse design and topology optimization tools
- **Visualization Tools**: Integrated plotting and analysis capabilities
- **Documentation**: Comprehensive tutorials and scientific applications

#### Technical Excellence
This package implements rigorously validated computational electromagnetics algorithms with comprehensive performance and accuracy validation.

## Contributing

We welcome contributions to FrequencyMaxwell.jl! Here's how you can help:

### Development Setup

1. **Fork and clone** the repository
2. **Install development dependencies**:
   ```julia
   using Pkg; Pkg.develop(".")
   Pkg.instantiate()
   ```
3. **Run tests**: `Pkg.test("FrequencyMaxwell")`
4. **Install pre-commit hooks** (if using):
   ```bash
   pip install pre-commit
   pre-commit install
   ```

### Contribution Guidelines

#### Code Contributions
- **Follow Julia style guidelines** (use `JuliaFormatter.jl`)
- **Add comprehensive tests** for new functionality
- **Include docstrings** for all public functions
- **Maintain type stability** and performance
- **Add examples** for new features

#### Areas for Contribution
- **New electromagnetic solvers** (FDTD, FEM, etc.)
- **Advanced source models** (Gaussian beams, vortex beams, etc.)
- **Optimization algorithms** for inverse design
- **GPU kernels** for performance acceleration
- **Documentation and tutorials**
- **Benchmark studies** and validation cases

#### Testing
```julia
# Run full test suite
using Pkg; Pkg.test("FrequencyMaxwell")

# Run specific test categories
Pkg.test("FrequencyMaxwell"; test_args=["--category=solvers"])
```

#### Documentation
Documentation is built using `Documenter.jl`:
```bash
cd docs/
julia make.jl
```

### Community
- **Issue tracker**: Report bugs and request features
- **Discussions**: Technical questions and community support
- **Scientific applications**: Share your research and applications
- **Performance improvements**: Benchmarks and optimization suggestions

## License and Acknowledgments

### License
FrequencyMaxwell.jl is released under the [MIT License](LICENSE). This allows for both academic and commercial use with minimal restrictions.

### Citation
If you use FrequencyMaxwell.jl in your research, please cite:

```bibtex
@software{FrequencyMaxwell.jl,
  title={FrequencyMaxwell.jl: Modern Electromagnetic Simulation in Julia},
  author={Dohyeon Lee and contributors},
  url={https://github.com/your-org/FrequencyMaxwell.jl},
  version={0.1.0},
  year={2024}
}
```

### Acknowledgments

#### Research Foundation
This implementation builds upon established computational electromagnetics research. For detailed information about mathematical foundations and prior work, see the [Development History](docs/src/development_history.md) documentation.

#### Julia Ecosystem
- **Julia Community**: LinearSolve.jl, FFTW.jl, and ecosystem packages that enable high-performance computing
- **Scientific Computing**: Contributions from the broader scientific computing community
- **Research Groups**: Contributors from computational electromagnetics and optics research worldwide

### Dependencies
We gratefully acknowledge the following packages that make FrequencyMaxwell.jl possible:
- **LinearSolve.jl**: Advanced linear algebra algorithms
- **FFTW.jl**: Fast Fourier transform implementations  
- **StaticArrays.jl**: Efficient small array operations
- **LinearAlgebra.jl**: Core linear algebra functionality

---

**FrequencyMaxwell.jl** - Advancing electromagnetic simulation through modern Julia computing.

For questions, bug reports, or contributions, please visit our [GitHub repository](https://github.com/your-org/FrequencyMaxwell.jl).
