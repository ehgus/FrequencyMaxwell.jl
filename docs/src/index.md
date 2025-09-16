# FrequencyMaxwell.jl

Modern Julia package for electromagnetic field simulation using frequency-domain Maxwell equations.

## Features

- **Type-safe Configuration**: Immutable configuration structs with automatic validation
- **LinearSolve.jl Integration**: Flexible linear algebra backends for performance
- **Memory Efficient**: Smart memory management with in-place operations
- **Extensible Architecture**: Modular design for easy extension and customization
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

## Future Implementation

FrequencyMaxwell.jl is actively developed with planned enhancements including GPU acceleration support and automatic differentiation capabilities for optimization workflows. These features are designed to maintain vendor-agnostic hardware support and seamless integration with the Julia ecosystem.
