# Phantoms and Geometry

FrequencyMaxwell.jl provides utilities for creating material distributions (phantoms) that represent the electromagnetic properties of objects to be simulated. These phantoms define the spatial variation of permittivity that creates electromagnetic scattering.

## Overview

In electromagnetic simulation, a phantom is a 3D array representing the spatial distribution of material properties. FrequencyMaxwell.jl uses relative permittivity phantoms where:

- **Background**: Regions with `permittivity_bg` (defined in configuration)
- **Objects**: Regions with different permittivity values creating scattering

The phantom array has the same dimensions as the computational grid defined in `ConvergentBornConfig`.

## Built-in Phantom Generators

### phantom_bead

Creates a spherical dielectric bead centered in the computational domain.

```julia
phantom = phantom_bead(grid_size, permittivity_values, radius)
```

**Parameters**:
- `grid_size`: `NTuple{3, Int}` - Grid dimensions (nx, ny, nz)
- `permittivity_values`: `Vector` - Permittivity values for each bead
- `radius`: `Real` - Radius in grid points

**Examples**:

```julia
# Single bead with n=1.5 (polystyrene-like)
config = ConvergentBornConfig(
    grid_size = (128, 128, 64),
    # ... other parameters
)

phantom_single = phantom_bead(
    config.grid_size,
    [1.5^2],  # Relative permittivity (n² = 2.25)
    16.0      # 16 grid points radius
)

# Multiple concentric beads (core-shell structure)
phantom_coreshell = phantom_bead(
    config.grid_size,
    [1.8^2, 1.4^2],  # Core: n=1.8, Shell: n=1.4
    [8.0, 16.0]      # Inner radius: 8, Outer radius: 16
)
```

### phantom_plate

Creates a flat plate or slab geometry.

```julia
phantom = phantom_plate(grid_size, permittivity_values, thickness, orientation)
```

**Parameters**:
- `grid_size`: `NTuple{3, Int}` - Grid dimensions
- `permittivity_values`: `Vector` - Permittivity values
- `thickness`: `Real` - Thickness in grid points
- `orientation`: `Symbol` - `:x`, `:y`, or `:z` for plate normal direction

**Examples**:

```julia
# Thin plate perpendicular to z-axis (xy-plane)
phantom_plate_z = phantom_plate(
    config.grid_size,
    [1.6^2],  # Glass-like material
    8.0,      # 8 grid points thick
    :z        # Normal along z-axis
)

# Thick slab perpendicular to x-axis
phantom_slab_x = phantom_plate(
    config.grid_size,
    [2.0^2],  # High index material
    32.0,     # 32 grid points thick
    :x        # Normal along x-axis
)
```

### phantom_cylinder

Creates a cylindrical geometry.

```julia
phantom = phantom_cylinder(grid_size, permittivity_values, radius, height, orientation)
```

**Parameters**:
- `grid_size`: `NTuple{3, Int}` - Grid dimensions
- `permittivity_values`: `Vector` - Permittivity values
- `radius`: `Real` - Cylinder radius in grid points
- `height`: `Real` - Cylinder height in grid points
- `orientation`: `Symbol` - `:x`, `:y`, or `:z` for cylinder axis direction

**Examples**:

```julia
# Cylinder along z-axis (fiber-like)
phantom_fiber = phantom_cylinder(
    config.grid_size,
    [1.45^2],  # Optical fiber core
    12.0,      # 12 grid points radius
    48.0,      # 48 grid points height
    :z         # Axis along z
)

# Short cylinder along x-axis
phantom_rod = phantom_cylinder(
    config.grid_size,
    [1.7^2],   # High index rod
    8.0,       # 8 grid points radius
    24.0,      # 24 grid points length
    :x         # Axis along x
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
config = ConvergentBornConfig(
    grid_size = (128, 128, 64),
    resolution = (50e-9, 50e-9, 50e-9),  # 50 nm voxels
    # ... other parameters
)

# Bead with 800 nm diameter (16 grid points radius)
radius_points = 16.0
radius_physical = radius_points * config.resolution[1]  # 800e-9 m

phantom = phantom_bead(config.grid_size, [1.5^2], radius_points)
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

phantom = phantom_bead(config.grid_size, [ε_complex], 16.0)
```

## Advanced Phantom Creation

### Custom Phantom Generation

Create custom phantoms using direct array manipulation:

```julia
function custom_phantom(grid_size, config)
    phantom = ones(ComplexF64, grid_size)  # Background permittivity = 1.0

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
            phantom[i, j, k] = 1.5^2 * (1 - 0.1 * r/20.0)
        end
    end

    return phantom
end

phantom = custom_phantom(config.grid_size, config)
```

### Combining Phantoms

Combine multiple phantom geometries:

```julia
# Create base phantom
phantom = ones(ComplexF64, config.grid_size)

# Add bead
bead = phantom_bead(config.grid_size, [1.5^2], 16.0)
phantom .= phantom .* (bead .== 1.0) + bead .* (bead .!= 1.0)

# Add plate
plate = phantom_plate(config.grid_size, [1.8^2], 4.0, :z)
phantom .= phantom .* (plate .== 1.0) + plate .* (plate .!= 1.0)
```

### Importing External Phantoms

Load phantoms from external sources:

```julia
using HDF5  # Or other file format libraries

# Load from HDF5 file
function load_phantom(filename, dataset_name)
    h5open(filename, "r") do file
        phantom = read(file, dataset_name)
        return Complex.(phantom)  # Ensure complex type
    end
end

phantom = load_phantom("my_phantom.h5", "permittivity")
```

## Validation and Best Practices

### Phantom Validation

Validate phantom properties before simulation:

```julia
function validate_phantom(phantom, config)
    # Check dimensions
    if size(phantom) != config.grid_size
        error("Phantom size mismatch")
    end

    # Check for reasonable permittivity values
    max_perm = maximum(real(phantom))
    min_perm = minimum(real(phantom))

    if max_perm > 20.0
        @warn "Very high permittivity detected: $max_perm"
    end

    if min_perm < 0.1
        @warn "Very low permittivity detected: $min_perm"
    end

    # Check for NaN or Inf values
    if any(!isfinite.(phantom))
        error("Phantom contains NaN or Inf values")
    end

    return true
end

validate_phantom(phantom, config)
```

### Performance Considerations

#### Memory Usage

Large phantoms consume significant memory:

```julia
# Float64 phantom memory usage
nx, ny, nz = config.grid_size
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
voxel_size = config.resolution[1]
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
background_perm = config.permittivity_bg
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
z_center = div(config.grid_size[3], 2)
phantom_xy = phantom[:, :, z_center]

# Extract central xz-plane
y_center = div(config.grid_size[2], 2)
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
    background = config.permittivity_bg
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
# Cell-like phantom (nucleus + cytoplasm)
function cell_phantom(grid_size, nucleus_permittivity, cytoplasm_permittivity,
                     cell_radius, nucleus_radius)
    phantom = ones(ComplexF64, grid_size)

    # Cell body
    cell = phantom_bead(grid_size, [cytoplasm_permittivity], cell_radius)

    # Nucleus
    nucleus = phantom_bead(grid_size, [nucleus_permittivity], nucleus_radius)

    # Combine
    phantom .= phantom .* (cell .== 1.0) + cell .* (cell .!= 1.0)
    phantom .= phantom .* (nucleus .== 1.0) + nucleus .* (nucleus .!= 1.0)

    return phantom
end

# Typical cell parameters
cell = cell_phantom(
    config.grid_size,
    1.38^2,  # Nucleus: n=1.38
    1.36^2,  # Cytoplasm: n=1.36
    20.0,    # Cell radius: 20 grid points
    8.0      # Nucleus radius: 8 grid points
)
```

### Optical Devices

```julia
# Waveguide structure
function waveguide_phantom(grid_size, core_permittivity, cladding_permittivity,
                          core_width, core_height)
    phantom = cladding_permittivity * ones(ComplexF64, grid_size)

    nx, ny, nz = grid_size
    cx, cy = div(nx, 2), div(ny, 2)

    # Create rectangular core
    x_range = (cx - div(core_width, 2)):(cx + div(core_width, 2))
    y_range = (cy - div(core_height, 2)):(cy + div(core_height, 2))

    phantom[x_range, y_range, :] .= core_permittivity

    return phantom
end

waveguide = waveguide_phantom(
    config.grid_size,
    1.45^2,  # Core: n=1.45
    1.44^2,  # Cladding: n=1.44
    8,       # Core width: 8 grid points
    8        # Core height: 8 grid points
)
```

These phantom utilities provide the foundation for electromagnetic simulations, enabling users to model complex geometries and material distributions with ease.