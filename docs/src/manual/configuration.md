# Solver Configuration

FrequencyMaxwell.jl uses an integrated configuration system where all parameters are directly embedded in the `ConvergentBornSolver` object. This streamlined approach provides type-safe parameter validation and automatic type promotion.

## ConvergentBornSolver Configuration

The `ConvergentBornSolver` contains all physical and numerical parameters needed for electromagnetic simulations.

### Basic Usage

```julia
using FrequencyMaxwell

solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64)
)
```

### Complete Parameter Reference

#### Physical Parameters

**`wavelength`** (required)
- **Type**: `T <: AbstractFloat`
- **Description**: Wavelength of electromagnetic radiation in vacuum
- **Units**: Meters
- **Validation**: Must be positive (`wavelength > 0`)
- **Example**: `500e-9` for 500 nm green light

**`permittivity_bg`** (required)
- **Type**: `T <: AbstractFloat`
- **Description**: Background relative permittivity (εᵣ = n²)
- **Units**: Dimensionless
- **Validation**: Must be positive (`permittivity_bg > 0`)
- **Example**: `1.33^2 ≈ 1.77` for water background

#### Spatial Discretization

**`resolution`** (required)
- **Type**: `NTuple{3, T}`
- **Description**: Voxel size in each spatial dimension
- **Units**: Meters
- **Validation**: All components must be positive
- **Example**: `(50e-9, 50e-9, 50e-9)` for 50 nm isotropic voxels

**`grid_size`** (required)
- **Type**: `NTuple{3, Int}`
- **Description**: Number of grid points in each spatial dimension
- **Units**: Dimensionless (grid points)
- **Validation**: All components must be positive
- **Example**: `(128, 128, 64)` for a 6.4 × 6.4 × 3.2 μm domain

#### Solver Parameters

**`tolerance`** (optional, default: `1e-6`)
- **Type**: `T <: AbstractFloat`
- **Description**: Convergence tolerance for iterative solver
- **Units**: Dimensionless
- **Validation**: Must be in range `[0, 1]`
- **Recommendation**: `1e-6` for most applications, `1e-8` for high precision

**`iterations_max`** (optional, default: `-1`)
- **Type**: `Int`
- **Description**: Maximum number of Born iterations (-1 for automatic)
- **Units**: Dimensionless
- **Validation**: Must be positive or -1
- **Recommendation**: Use `-1` for automatic determination, or set a specific limit

**`linear_solver`** (optional, default: `KrylovJL_GMRES()`)
- **Type**: `SciMLLinearSolveAlgorithm`
- **Description**: LinearSolve.jl algorithm for linear system solution
- **Options**: Any LinearSolve.jl algorithm (KrylovJL_GMRES, BiCGStabL, etc.)
- **Example**: `KrylovJL_GMRES()` for recommended default choice

#### Hardware Configuration

**`device`** (optional, default: `:cpu`)
- **Type**: `Symbol`
- **Description**: Device backend for computation
- **Options**: `:cpu` (multi-threaded), `:cuda` (NVIDIA), `:amdgpu` (AMD), `:metal` (Apple), `:oneapi` (Intel)
- **Validation**: Requires corresponding GPU package to be loaded (CUDA.jl, AMDGPU.jl, Metal.jl, or oneAPI.jl)
- **Example**: `device = :cuda` for NVIDIA GPU acceleration

#### Boundary Configuration

**`boundary_thickness`** (optional, default: `(0.0, 0.0, 0.0)`)
- **Type**: `NTuple{3, T}`
- **Description**: PML boundary layer thickness in each direction
- **Units**: Meters
- **Usage**: Set to non-zero values for absorbing boundary conditions

**`field_attenuation`** (optional, default: `(0.0, 0.0, 0.0)`)
- **Type**: `NTuple{3, T}`
- **Description**: Field attenuation layer thickness in each direction
- **Units**: Meters

**`field_attenuation_sharpness`** (optional, default: `1.0`)
- **Type**: `T <: AbstractFloat`
- **Description**: Sharpness factor for field attenuation (0-1)
- **Units**: Dimensionless

**`periodic_boundary`** (optional, default: `(true, true, false)`)
- **Type**: `NTuple{3, Bool}`
- **Description**: Periodic boundary conditions (x, y, z)
- **Recommendation**: `(true, true, false)` for typical 3D simulations


## Configuration Utilities

### Helper Functions

```julia
# Get total domain size
domain = domain_size(solver)  # Returns (Lx, Ly, Lz) in meters

# Get voxel spacing (same as resolution)
spacing = grid_spacing(solver)  # Returns (dx, dy, dz) in meters

# Get background wavenumber
k0 = wavenumber_background(solver)  # Returns k₀ = 2π/λ in rad/m
```

### Validation and Type Promotion

The configuration system automatically:

1. **Validates all parameters** with descriptive error messages
2. **Promotes types** to ensure consistency (e.g., Int to Float64)
3. **Provides meaningful defaults** for optional parameters

```julia
# This automatically promotes Int to Float64
config = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 2,        # Int promoted to Float64
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64)
)

# Type information is preserved
println(typeof(config))  # ConvergentBornConfig{Float64}
```

## Advanced Configuration

### GPU-Accelerated Simulations

For large-scale problems, GPU acceleration provides significant speedup:

```julia
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (256, 256, 128),
    device = :cuda,              # Use NVIDIA GPU
    tolerance = 1e-6
)

# For AMD GPUs: device = :amdgpu
# For Apple Silicon: device = :metal
# For Intel GPUs: device = :oneapi
```

### High-Resolution Simulations

For high-resolution simulations with large grids:

```julia
config = ConvergentBornConfig{Float32}(  # Memory-efficient precision
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (25e-9, 25e-9, 25e-9),  # Higher resolution
    grid_size = (256, 256, 128)          # Larger grids
)
```

### Large Simulations

For large simulations, consider memory-efficient configurations:

```julia
config = ConvergentBornConfig{Float32}(  # Use Float32 for memory efficiency
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (512, 512, 256),
    tolerance = 1e-5               # Slightly relaxed tolerance for speed
)
```

### High-Precision Problems

For applications requiring high numerical accuracy:

```julia
config = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64),
    tolerance = 1e-7,              # Higher precision
    linear_solver = KrylovJL_GMRES()  # Recommended default choice
)
```

### Custom Linear Solvers

FrequencyMaxwell.jl supports any LinearSolve.jl algorithm:

```julia
using LinearSolve

# Different solver algorithms
config_gmres = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64),
    linear_solver = KrylovJL_GMRES(),
    tolerance = 1e-6
)

config_bicgstab = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64),
    linear_solver = KrylovJL_GMRES()  # Recommended default
    tolerance = 1e-6
)

config_idr = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64),
    linear_solver = KrylovJL_GMRES()  # Recommended default
    tolerance = 1e-6
)
```

## Best Practices

### Development Workflow

1. **Start small**: Begin with small grid sizes (`(64, 64, 32)`) for development
2. **Check convergence**: Monitor solver convergence with different parameters
3. **Scale up gradually**: Increase resolution and grid size for production runs

### Parameter Selection Guidelines

- **Wavelength**: Use the vacuum wavelength of your electromagnetic source
- **Resolution**: Choose 10-20 points per wavelength in the background medium
- **Grid size**: Ensure objects of interest are well-sampled and boundaries are far enough
- **Tolerance**: `1e-6` for most applications, `1e-8` for high-precision work

### Performance Optimization

- **Use appropriate precision**: Consider Float32 for memory-limited problems
- **Algorithm choice**: Try different linear solvers if convergence is slow
- **Boundary configuration**: Optimize boundary thickness and field attenuation for your specific problem

## Configuration Display

Configurations have custom `show` methods for clear REPL display:

```julia
julia> config = ConvergentBornConfig(wavelength=500e-9, permittivity_bg=1.77,
                                    resolution=(50e-9, 50e-9, 50e-9), grid_size=(128, 128, 64))

ConvergentBornConfig{Float64}:
  Physical parameters:
    wavelength: 5.0e-7 m
    permittivity_bg: 1.77
    resolution: (5.0e-8, 5.0e-8, 5.0e-8) m
    grid_size: (128, 128, 64)
  Domain: 6.4 × 6.4 × 3.2 μm
  Solver: KrylovJL_GMRES(), tolerance=1.0e-6, iterations_max=-1
```

This provides immediate insight into the physical and numerical parameters of your simulation.