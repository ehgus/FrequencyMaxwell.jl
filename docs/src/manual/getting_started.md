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

1. **Configure the solver** with physical and numerical parameters
2. **Create a solver instance** from the configuration
3. **Define electromagnetic sources** (plane waves, multiple plane waves)
4. **Generate or load a material distribution** (permittivity phantom)
5. **Solve the electromagnetic scattering problem**
6. **Analyze the resulting fields**

## Your First Simulation

Let's start with a simple scattering simulation of a dielectric bead:

```julia
using FrequencyMaxwell

# Step 1: Configure the solver
config = ConvergentBornConfig(
    wavelength = 500e-9,          # 500 nm wavelength (green light)
    permittivity_bg = 1.33^2,     # Water background (n = 1.33)
    resolution = (50e-9, 50e-9, 50e-9),  # 50 nm voxel size
    grid_size = (128, 128, 32),   # 6.4 × 6.4 × 1.6 μm computational domain
    tolerance = 1e-6              # Convergence tolerance
)

# Step 2: Create solver instance
solver = ConvergentBornSolver(config)

# Step 3: Define electromagnetic source
source = PlaneWaveSource(
    wavelength = config.wavelength,
    polarization = [1.0, 0.0, 0.0],  # x-polarized light
    k_vector = [0.0, 0.0, 1.0],      # propagating in +z direction
    amplitude = 1.0,                  # unit amplitude
    phase = 0.0                       # zero phase
)

# Step 4: Generate material distribution (dielectric bead)
permittivity_phantom = phantom_bead(
    config.grid_size,     # grid dimensions
    [1.5^2],              # permittivity of bead (n = 1.5)
    16.0                  # radius in grid points
)

# Step 5: Solve the electromagnetic problem
E_field, H_field = solve(solver, source, permittivity_phantom)

# Step 6: Analyze results
println("Electric field dimensions: ", size(E_field.E))
println("Magnetic field dimensions: ", size(H_field.H))

# Calculate field energy
total_energy = field_energy(E_field) + field_energy(H_field)
println("Total field energy: ", total_energy)
```

## Understanding the Results

The `solve` function returns two `ElectromagneticField` objects:

- `E_field`: Contains the total electric field (incident + scattered)
- `H_field`: Contains the total magnetic field (incident + scattered)

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
- **Memory Management**: Efficient memory usage for large-scale simulations
- **Precision**: Supports both Float32 and Float64 precision

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