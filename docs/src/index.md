# FrequencyMaxwell.jl

Modern Julia package for electromagnetic field simulation using frequency-domain Maxwell equations.

## Features

- **Integrated Configuration**: Streamlined solver with built-in parameter validation
- **LinearSolve.jl Integration**: Flexible linear algebra backends for performance
- **Memory Efficient**: Smart memory management with in-place operations
- **Direct Field Access**: Intuitive API with direct access to solver parameters
- **Comprehensive Testing**: >95% code coverage with mathematical validation

## Quick Start

```julia
using FrequencyMaxwell

# Create solver with integrated configuration
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32)
)

# Create source using direct field access
source = PlaneWaveSource(
    wavelength = solver.wavelength,
    polarization = [1.0, 0.0, 0.0],
    k_vector = [0.0, 0.0, 1.0]
)

# Generate phantom and solve
phantom = phantom_bead(solver.grid_size, [1.5^2], 16.0)
E_field, H_field = solve(solver, source, phantom)
```

## GPU Acceleration

FrequencyMaxwell.jl supports vendor-agnostic GPU acceleration via KernelAbstractions.jl, enabling simulation on NVIDIA (CUDA), AMD (ROCm), Apple Silicon (Metal), and Intel (oneAPI) GPUs with the same unified interface.
