# Phantoms and Geometry

FrequencyMaxwell.jl provides utilities for creating material distributions (phantoms) that represent the electromagnetic properties of objects to be simulated. These phantoms define the spatial variation of permittivity that creates electromagnetic scattering.

## Overview

In electromagnetic simulation, a phantom represents the spatial distribution of material properties. FrequencyMaxwell.jl uses the `Medium` type to encapsulate both the permittivity distribution and background permittivity:

- **Medium object**: Contains both the 3D permittivity array and the background permittivity value
- **Background**: Regions with `permittivity_bg`
- **Objects**: Regions with different permittivity values creating scattering

The phantom's permittivity array has the same dimensions as the computational grid defined in `ConvergentBornSolver`. All built-in phantom generators return `Medium` objects ready for use with the solver.

## Built-in Phantom Generators

### phantom_bead

Creates a spherical dielectric bead centered in the computational domain, returning a `Medium` object.

```julia
medium = phantom_bead(grid_size, permittivity_values, radius; permittivity_bg=1.0)
```

**Parameters**:
- `grid_size`: `NTuple{3, Int}` - Grid dimensions (nx, ny, nz)
- `permittivity_values`: `Vector` - Permittivity values for each bead
- `radius`: `Real` - Radius in grid points
- `permittivity_bg`: `Real` - Background permittivity (keyword argument)

**Returns**: `Medium` object containing the permittivity distribution and background

**Examples**:

```julia
# Single bead with n=1.5 (polystyrene-like) in water
solver = ConvergentBornSolver(
    grid_size = (128, 128, 64),
    # ... other parameters
)

medium = phantom_bead(
    solver.grid_size,
    [1.5^2],              # Relative permittivity (n² = 2.25)
    16.0,                 # 16 grid points radius
    permittivity_bg = 1.33^2  # Water background
)

# Use directly with solver
EMfield = solve(solver, source, medium)

# Multiple concentric beads (core-shell structure)
medium_coreshell = phantom_bead(
    solver.grid_size,
    [1.8^2, 1.4^2],      # Core: n=1.8, Shell: n=1.4
    [8.0, 16.0],         # Inner radius: 8, Outer radius: 16
    permittivity_bg = 1.0
)
```

### phantom_plate

Creates a flat plate or slab geometry, returning a `Medium` object.

```julia
medium = phantom_plate(grid_size, permittivity_values, thickness; permittivity_bg=1.0)
```

**Parameters**:
- `grid_size`: `NTuple{3, Int}` - Grid dimensions
- `permittivity_values`: `Vector` - Permittivity values
- `thickness`: `Real` - Thickness in grid points
- `permittivity_bg`: `Real` - Background permittivity (keyword argument)

**Returns**: `Medium` object containing the permittivity distribution and background

**Examples**:

```julia
# Thin plate in water
medium_plate = phantom_plate(
    solver.grid_size,
    [1.6^2],              # Glass-like material
    8.0,                  # 8 grid points thick
    permittivity_bg = 1.33^2
)

# Thick slab in vacuum
medium_slab = phantom_plate(
    solver.grid_size,
    [2.0^2],              # High index material
    32.0,                 # 32 grid points thick
    permittivity_bg = 1.0
)
```

## Phantom Properties and Usage

### Coordinate System

Phantoms use a grid-based coordinate system where:
- Origin is at the center of the computational domain
- Grid points are indexed from 1 to grid_size[i] in each dimension
- Coordinates can be converted to physical units using the resolution

### Physical Size Calculations

```julia
solver = ConvergentBornSolver(
    grid_size = (128, 128, 64),
    resolution = (50e-9, 50e-9, 50e-9),  # 50 nm voxels
    # ... other parameters
)

# Bead with 800 nm diameter (16 grid points radius)
radius_points = 16.0
radius_physical = radius_points * solver.resolution[1]  # 800e-9 m

medium = phantom_bead(
    solver.grid_size,
    [1.5^2],
    radius_points,
    permittivity_bg = 1.33^2
)
```

### Material Properties

#### Relative Permittivity

Phantoms store relative permittivity (εᵣ), related to refractive index by:

```julia
n = 1.5  # Refractive index
εᵣ = n^2  # Relative permittivity = 2.25
```

Common materials:
- **Air/Vacuum**: εᵣ = 1.0 (n = 1.0)
- **Water**: εᵣ = 1.77 (n = 1.33)
- **Glass**: εᵣ = 2.25-2.56 (n = 1.5-1.6)
- **Polystyrene**: εᵣ = 2.56 (n = 1.6)
- **Silicon**: εᵣ ≈ 12 (n ≈ 3.5, wavelength dependent)

#### Loss and Absorption

For lossy materials, use complex permittivity:

```julia
# Material with loss (imaginary part)
n_real = 1.5
n_imag = 0.1  # Absorption
ε_complex = (n_real + 1im * n_imag)^2

medium = phantom_bead(
    solver.grid_size,
    [ε_complex],
    16.0,
    permittivity_bg = 1.33^2
)
```

## Advanced Phantom Creation

### Custom Phantom Generation

Create custom phantoms using direct array manipulation and wrap them in a `Medium`:

```julia
function custom_phantom(grid_size, permittivity_bg)
    # Create permittivity array
    perm_array = fill(Complex{Float64}(permittivity_bg), grid_size)

    # Get grid coordinates
    nx, ny, nz = grid_size
    x = range(-nx/2, nx/2, length=nx)
    y = range(-ny/2, ny/2, length=ny)
    z = range(-nz/2, nz/2, length=nz)

    # Create custom structure
    for (i, xi) in enumerate(x), (j, yj) in enumerate(y), (k, zk) in enumerate(z)
        r = sqrt(xi^2 + yj^2 + zk^2)

        # Example: Graded index structure
        if r < 20.0
            perm_array[i, j, k] = 1.5^2 * (1 - 0.1 * r/20.0)
        end
    end

    # Wrap in Medium object
    return Medium(perm_array, permittivity_bg)
end

medium = custom_phantom(solver.grid_size, 1.33^2)
EMfield = solve(solver, source, medium)
```

### Combining Phantoms

Combine multiple phantom geometries by extracting and combining their permittivity arrays:

```julia
using FrequencyMaxwell: permittivity

# Start with background
permittivity_bg = 1.33^2
perm_combined = fill(Complex{Float64}(permittivity_bg), solver.grid_size)

# Add bead
medium_bead = phantom_bead(solver.grid_size, [1.5^2], 16.0, permittivity_bg=permittivity_bg)
perm_bead = permittivity(medium_bead)
mask_bead = abs.(perm_bead .- permittivity_bg) .> 1e-10
perm_combined[mask_bead] .= perm_bead[mask_bead]

# Add plate
medium_plate = phantom_plate(solver.grid_size, [1.8^2], 4.0, permittivity_bg=permittivity_bg)
perm_plate = permittivity(medium_plate)
mask_plate = abs.(perm_plate .- permittivity_bg) .> 1e-10
perm_combined[mask_plate] .= perm_plate[mask_plate]

# Create final Medium
medium_combined = Medium(perm_combined, permittivity_bg)
```

### Importing External Phantoms

Load phantoms from external sources and wrap them in Medium objects:

```julia
using HDF5  # Or other file format libraries

# Load from HDF5 file
function load_phantom_medium(filename, dataset_name, permittivity_bg)
    perm_array = h5open(filename, "r") do file
        data = read(file, dataset_name)
        return Complex{Float64}.(data)  # Ensure complex type
    end

    return Medium(perm_array, permittivity_bg)
end

medium = load_phantom_medium("my_phantom.h5", "permittivity", 1.33^2)
EMfield = solve(solver, source, medium)
```

## Validation and Best Practices

### Phantom Validation

Validate Medium properties before simulation:

```julia
using FrequencyMaxwell: permittivity, grid_size

function validate_medium(medium, solver)
    # Check dimensions
    if grid_size(medium) != solver.grid_size
        error("Medium size mismatch: got $(grid_size(medium)), expected $(solver.grid_size)")
    end

    # Get permittivity array
    perm = permittivity(medium)

    # Check for reasonable permittivity values
    max_perm = maximum(real(perm))
    min_perm = minimum(real(perm))

    if max_perm > 20.0
        @warn "Very high permittivity detected: $max_perm"
    end

    if min_perm < 0.1
        @warn "Very low permittivity detected: $min_perm"
    end

    # Check for NaN or Inf values
    if any(!isfinite.(perm))
        error("Medium contains NaN or Inf values")
    end

    return true
end

validate_medium(medium, solver)
```

### Performance Considerations

#### Memory Usage

Large phantoms consume significant memory:

```julia
# Float64 phantom memory usage
nx, ny, nz = solver.grid_size
memory_gb = nx * ny * nz * 16 / 1e9  # 16 bytes per Complex{Float64}
println("Phantom memory: $(memory_gb) GB")

# Use Float32 for memory savings
phantom_f32 = ComplexF32.(phantom)  # 50% memory reduction
```

#### Grid Resolution

Choose appropriate resolution for phantom features:

```julia
# Ensure features are adequately sampled
feature_size_physical = 200e-9  # 200 nm feature
voxel_size = solver.resolution[1]
points_per_feature = feature_size_physical / voxel_size

if points_per_feature < 4
    @warn "Feature may be undersampled: $points_per_feature points per feature"
end
```

### Physical Realism

#### Contrast Limitations

Very high contrast phantoms may cause convergence issues:

```julia
# Calculate contrast
background_perm = solver.permittivity_bg
max_contrast = maximum(real(phantom)) / background_perm

if max_contrast > 4.0
    @warn "High contrast detected: $max_contrast. Consider using smaller tolerance."
end
```

#### Smooth Transitions

Smooth permittivity transitions often improve convergence:

```julia
function smooth_phantom(phantom, σ=1.0)
    # Apply Gaussian smoothing (requires Images.jl or similar)
    using ImageFiltering
    return imfilter(phantom, Kernel.gaussian(σ))
end

phantom_smooth = smooth_phantom(phantom, 0.5)
```

## Phantom Analysis and Visualization

### Cross-sections

Extract 2D cross-sections for visualization:

```julia
# Extract central xy-plane
z_center = div(solver.grid_size[3], 2)
phantom_xy = phantom[:, :, z_center]

# Extract central xz-plane
y_center = div(solver.grid_size[2], 2)
phantom_xz = phantom[:, y_center, :]
```

### Statistics

Analyze phantom properties:

```julia
function phantom_statistics(phantom, config)
    # Basic statistics
    mean_perm = mean(real(phantom))
    std_perm = std(real(phantom))

    # Contrast statistics
    background = solver.permittivity_bg
    contrast = real(phantom) ./ background
    max_contrast = maximum(contrast)
    min_contrast = minimum(contrast)

    println("Mean permittivity: $mean_perm")
    println("Standard deviation: $std_perm")
    println("Contrast range: $min_contrast to $max_contrast")

    # Volume fractions
    high_index = count(real(phantom) .> background * 1.1)
    total_voxels = length(phantom)
    volume_fraction = high_index / total_voxels

    println("High index volume fraction: $(volume_fraction * 100)%")
end

phantom_statistics(phantom, config)
```

## Common Phantom Types

### Biological Cells

```julia
using FrequencyMaxwell: permittivity

# Cell-like phantom (nucleus + cytoplasm)
function cell_phantom_medium(grid_size, nucleus_permittivity, cytoplasm_permittivity,
                            cell_radius, nucleus_radius, permittivity_bg)
    # Start with background
    perm_array = fill(Complex{Float64}(permittivity_bg), grid_size)

    # Add cell body
    cell = phantom_bead(grid_size, [cytoplasm_permittivity], cell_radius,
                       permittivity_bg=permittivity_bg)
    perm_cell = permittivity(cell)
    mask_cell = abs.(perm_cell .- permittivity_bg) .> 1e-10
    perm_array[mask_cell] .= perm_cell[mask_cell]

    # Add nucleus (overwrites cell center)
    nucleus = phantom_bead(grid_size, [nucleus_permittivity], nucleus_radius,
                          permittivity_bg=permittivity_bg)
    perm_nucleus = permittivity(nucleus)
    mask_nucleus = abs.(perm_nucleus .- permittivity_bg) .> 1e-10
    perm_array[mask_nucleus] .= perm_nucleus[mask_nucleus]

    return Medium(perm_array, permittivity_bg)
end

# Typical cell parameters
medium_cell = cell_phantom_medium(
    solver.grid_size,
    1.38^2,  # Nucleus: n=1.38
    1.36^2,  # Cytoplasm: n=1.36
    20.0,    # Cell radius: 20 grid points
    8.0,     # Nucleus radius: 8 grid points
    1.33^2   # Water background
)
```

### Optical Devices

```julia
# Waveguide structure
function waveguide_phantom_medium(grid_size, core_permittivity, cladding_permittivity,
                                 core_width, core_height)
    # Start with cladding
    perm_array = fill(Complex{Float64}(cladding_permittivity), grid_size)

    nx, ny, nz = grid_size
    cx, cy = div(nx, 2), div(ny, 2)

    # Create rectangular core
    x_range = (cx - div(core_width, 2)):(cx + div(core_width, 2))
    y_range = (cy - div(core_height, 2)):(cy + div(core_height, 2))

    perm_array[x_range, y_range, :] .= core_permittivity

    # Use cladding as background
    return Medium(perm_array, cladding_permittivity)
end

medium_waveguide = waveguide_phantom_medium(
    solver.grid_size,
    1.45^2,  # Core: n=1.45
    1.44^2,  # Cladding: n=1.44
    8,       # Core width: 8 grid points
    8        # Core height: 8 grid points
)
```

These phantom utilities provide the foundation for electromagnetic simulations, enabling users to model complex geometries and material distributions with ease. All phantom functions return `Medium` objects that are directly compatible with the solver's `solve()` function.