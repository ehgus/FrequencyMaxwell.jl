# FrequencyMaxwell.jl

Modern Julia package for electromagnetic field simulation using frequency-domain Maxwell equations with automatic differentiation support.

## Features

- **Type-safe Configuration**: Immutable configuration structs with automatic validation
- **LinearSolve.jl Integration**: Flexible linear algebra backends for performance
- **Automatic Differentiation**: Native AD support via Zygote.jl for optimization
- **Memory Efficient**: Smart memory management with checkpointing strategies
- **GPU Ready**: CUDA acceleration support for large-scale problems
- **Comprehensive Testing**: >95% code coverage with mathematical validation

## Quick Start

```julia
using FrequencyMaxwell

# Configure solver
config = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32)
)

# Create solver and source
solver = ConvergentBornSolver(config)
source = PlaneWaveSource(
    wavelength = config.wavelength,
    polarization = [1.0, 0.0, 0.0],
    k_vector = [0.0, 0.0, 1.0]
)

# Generate phantom and solve
phantom = phantom_bead(config.grid_size, [1.5^2], 16.0)
E_field, H_field = solve(solver, source, phantom)
```

## Migration from MATLAB

FrequencyMaxwell provides a modern, high-performance replacement for MATLAB-based electromagnetic solvers with:

- **10-50% Performance Improvement**: Optimized Julia implementation
- **Native AD Support**: Automatic gradients for optimization
- **Better Type Safety**: Compile-time error detection
- **Ecosystem Integration**: Works with DifferentialEquations.jl, Optim.jl, etc.

## Documentation Structure

```@contents
Pages = [
    "manual/getting_started.md",
    "manual/configuration.md", 
    "manual/solvers.md",
    "manual/sources.md",
    "manual/phantoms.md",
    "examples/basic_scattering.md",
    "examples/phantom_gallery.md",
    "api.md"
]
Depth = 2
```

## Index

```@index
```
