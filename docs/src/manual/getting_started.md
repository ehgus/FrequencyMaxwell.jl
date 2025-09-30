# Getting Started

This guide will help you get up and running with FrequencyMaxwell.jl for electromagnetic simulation.

## Installation

FrequencyMaxwell.jl can be installed using package manager:

```julia
using Pkg
Pkg.add("https://github.com/ehgus/FrequencyMaxwell.jl.git")
```

## Basic Workflow

The typical workflow for electromagnetic simulation with FrequencyMaxwell.jl follows these steps:

1. **Create a solver** with integrated configuration
2. **Define electromagnetic sources** (plane waves, multiple plane waves)
3. **Generate or load a material distribution** (permittivity phantom)
4. **Solve the electromagnetic scattering problem**
5. **Analyze the resulting fields**

## Your First Simulation

Let's start with a simple scattering simulation of a dielectric bead:

```julia
using FrequencyMaxwell

# Step 1: Create solver with integrated configuration
solver = ConvergentBornSolver(
    wavelength = 500e-9,          # 500 nm wavelength (green light)
    permittivity_bg = 1.33^2,     # Water background (n = 1.33)
    resolution = (50e-9, 50e-9, 50e-9),  # 50 nm voxel size
    grid_size = (128, 128, 32),   # 6.4 × 6.4 × 1.6 μm computational domain
    tolerance = 1e-6              # Convergence tolerance
)

# Step 2: Define electromagnetic source
source = PlaneWaveSource(
    wavelength = solver.wavelength,
    polarization = [1.0, 0.0, 0.0],  # x-polarized light
    k_vector = [0.0, 0.0, 1.0],      # propagating in +z direction
    amplitude = 1.0,                  # unit amplitude
    phase = 0.0                       # zero phase
)

# Step 3: Generate material distribution (dielectric bead)
permittivity_phantom = phantom_bead(
    solver.grid_size,     # grid dimensions
    [1.5^2],              # permittivity of bead (n = 1.5)
    16.0                  # radius in grid points
)

# Step 4: Solve the electromagnetic problem
Efield, Hfield = solve(solver, source, permittivity_phantom)

# Step 5: Analyze results
println("Electric field dimensions: ", size(Efield))
println("Magnetic field dimensions: ", size(Hfield))

# Calculate field energy
total_energy = field_energy(Efield) + field_energy(Hfield)
println("Total field energy: ", total_energy)
```

## Understanding the Results

The `solve` function returns two `ElectromagneticField` objects:

- `Efield`: Contains the total electric field (incident + scattered)
- `Hfield`: Contains the total magnetic field (incident + scattered)

Each field object has:
- `.E` or `.H`: The 4D field array with dimensions `(3, nx, ny, nz)` representing the vector field components
- Methods for field analysis like `field_energy()`, `field_intensity()`, and `poynting_vector()`

## Key Concepts

### Configuration-First Design

FrequencyMaxwell.jl uses an immutable configuration approach where all simulation parameters are specified upfront in the `ConvergentBornConfig`. This provides:

- **Type safety**: Automatic validation of all parameters
- **Reproducibility**: Complete parameter specifications
- **Performance**: Compile-time optimizations based on configuration

### Hardware Flexibility

The package is designed to work across different hardware:

- **CPU**: Standard Julia arrays for development and small problems
- **GPU**: Vendor-agnostic GPU acceleration (NVIDIA CUDA, AMD ROCm, Apple Metal, Intel oneAPI)
- **Memory Management**: Efficient memory usage for large-scale simulations
- **Precision**: Supports both Float32 and Float64 precision

### GPU Acceleration

Enable GPU acceleration by specifying the `device` parameter:

```julia
# Using NVIDIA GPU (requires CUDA.jl)
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32),
    device = :cuda
)

# Other supported devices: :amdgpu (ROCm), :metal (Apple), :oneapi (Intel)
```

The package will automatically check for GPU availability and provide installation instructions if needed. Use `:cpu` for CPU-only computation.

### LinearSolve.jl Integration

FrequencyMaxwell.jl leverages the LinearSolve.jl ecosystem for maximum flexibility in linear algebra:

```julia
# Example: Using a different linear solver algorithm
config = ConvergentBornConfig(
    # ... other parameters ...
    linear_solver = KrylovJL_GMRES()  # Recommended default algorithm
)
```

## Common Gotchas

### Memory Usage

Large 3D simulations can use significant memory. Monitor memory usage and consider:
- Reducing grid size for development
- Using Float32 precision for memory savings
- Using memory-efficient precision for larger problems

### Convergence Issues

If the solver doesn't converge:
- Try adjusting `iterations_max` parameter (-1 for automatic, or set a specific limit)
- Adjust `tolerance` (increase for faster convergence, decrease for more precision)
- Try different `linear_solver` algorithms (KrylovJL_GMRES, BiCGStabL, etc.)
- Check that your permittivity phantom is reasonable

### Physical Units

FrequencyMaxwell.jl uses SI units throughout:
- Wavelength in meters (e.g., `500e-9` for 500 nm)
- Resolution in meters (e.g., `50e-9` for 50 nm)
- Permittivity as relative permittivity (dimensionless)

The package handles unit consistency automatically, but ensure your inputs use SI units.