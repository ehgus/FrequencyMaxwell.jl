# Phantom Gallery Example

This example showcases the phantom generation capabilities of FrequencyMaxwell.jl, demonstrating how to create various material distributions for electromagnetic simulations. The gallery covers common geometries and materials encountered in optical microscopy, flow cytometry, and materials characterization.

## Overview

Phantoms in FrequencyMaxwell.jl represent 3D distributions of material properties (relative permittivity) that create electromagnetic scattering. This example demonstrates:

- **Basic geometric shapes**: Spheres, plates, cylinders
- **Material properties**: Realistic optical materials
- **Multi-object phantoms**: Complex experimental scenarios
- **Parameter analysis**: Statistical characterization of phantoms

## Complete Example Code

```julia
using FrequencyMaxwell

function phantom_gallery_example()
    println("FrequencyMaxwell Phantom Gallery")
    println("=" ^ 35)

    # Common parameters for all phantoms
    grid_size = (128, 128, 64)
    resolution = (50e-9, 50e-9, 100e-9)  # 50×50×100 nm voxels

    println("Grid parameters:")
    println("  Grid size: $grid_size")
    println("  Resolution: $(resolution .* 1e9) nm")
    println("  Domain size: $((grid_size .* resolution) .* 1e6) μm")
    println()

    # 1. Single polystyrene bead
    println("1. Single Polystyrene Bead")
    println("-" ^ 28)

    bead_single = phantom_bead(
        grid_size,
        [1.59^2],        # Polystyrene n=1.59, ε=2.53
        12.0             # 12-pixel radius = 600 nm
    )

    analyze_phantom(bead_single, "Single bead", grid_size, resolution)
    println()

    # 2. Multiple beads with different materials
    println("2. Multiple Material Beads")
    println("-" ^ 27)

    # Create multiple beads with different materials
    beads = []

    # Silica bead
    bead_silica = phantom_bead(grid_size, [1.46^2], 8.0)
    push!(beads, ("Silica bead", bead_silica))

    # Polystyrene bead
    bead_ps = phantom_bead(grid_size, [1.59^2], 10.0)
    push!(beads, ("Polystyrene bead", bead_ps))

    # High-index bead (TiO2)
    bead_tio2 = phantom_bead(grid_size, [2.4^2], 6.0)
    push!(beads, ("TiO2 bead", bead_tio2))

    for (name, bead) in beads
        analyze_phantom(bead, name, grid_size, resolution)
    end
    println()

    # 3. Plate geometries
    println("3. Plate Geometries")
    println("-" ^ 19)

    # Glass cover slip (horizontal)
    plate_glass = phantom_plate(
        grid_size,
        [1.52^2],        # Borosilicate glass n=1.52
        16.0,            # 16-pixel thickness = 1.6 μm
        :z               # Normal along z-axis (horizontal plate)
    )

    # Thin film (vertical)
    film_thin = phantom_plate(
        grid_size,
        [2.4^2],         # High-index coating n=2.4
        4.0,             # 4-pixel thickness = 200 nm
        :y               # Normal along y-axis (vertical film)
    )

    analyze_phantom(plate_glass, "Glass cover slip", grid_size, resolution)
    analyze_phantom(film_thin, "Thin film coating", grid_size, resolution)
    println()

    # 4. Cylindrical geometries
    println("4. Cylindrical Geometries")
    println("-" ^ 24)

    # Optical fiber
    fiber = phantom_cylinder(
        grid_size,
        [1.46^2],        # Silica core n=1.46
        8.0,             # 8-pixel radius = 400 nm
        40.0,            # 40-pixel height = 4 μm
        :z               # Axis along z-direction
    )

    # Biological filament
    filament = phantom_cylinder(
        grid_size,
        [1.37^2],        # Protein-like n=1.37
        3.0,             # 3-pixel radius = 150 nm
        50.0,            # 50-pixel height = 5 μm
        :z               # Axis along z-direction
    )

    analyze_phantom(fiber, "Optical fiber", grid_size, resolution)
    analyze_phantom(filament, "Biological filament", grid_size, resolution)
    println()

    # 5. Material properties reference
    println("5. Material Properties Reference")
    println("-" ^ 33)

    materials = [
        ("Air/Vacuum", 1.0^2),
        ("Water", 1.33^2),
        ("Cytoplasm", 1.37^2),
        ("Cell nucleus", 1.39^2),
        ("Silica (SiO₂)", 1.46^2),
        ("Borosilicate glass", 1.52^2),
        ("Polystyrene", 1.59^2),
        ("Silicon (visible)", 4.0^2),
        ("Titanium dioxide", 2.4^2)
    ]

    println("Optical materials database:")
    println("Material                 n      ε      Δε (vs water)")
    println("-" ^ 50)

    water_perm = 1.33^2
    for (name, perm) in materials
        n = sqrt(real(perm))
        contrast = real(perm) / water_perm - 1
        println("$(rpad(name, 20)) $(rpad(round(n, digits=3), 6)) $(rpad(round(real(perm), digits=3), 6)) $(round(contrast, digits=3))")
    end
    println()

    # 6. Complex multi-object phantom
    println("6. Complex Multi-Object Phantom")
    println("-" ^ 33)

    complex_phantom = create_complex_phantom(grid_size)
    analyze_phantom(complex_phantom, "Complex scene", grid_size, resolution)
    println()

    # 7. Biological phantoms
    println("7. Biological Cell Phantoms")
    println("-" ^ 29)

    # Simple cell model
    cell_simple = create_cell_phantom(grid_size, :simple)
    analyze_phantom(cell_simple, "Simple cell", grid_size, resolution)

    # Cell with organelles
    cell_complex = create_cell_phantom(grid_size, :complex)
    analyze_phantom(cell_complex, "Cell with organelles", grid_size, resolution)
    println()

    println("Gallery Summary:")
    println("✓ All phantoms use realistic optical material properties")
    println("✓ Geometries suitable for electromagnetic scattering simulations")
    println("✓ Applications demonstrated:")
    println("  • Particle characterization (beads)")
    println("  • Optical component modeling (plates, fibers)")
    println("  • Biological imaging (cells, organelles)")
    println("  • Materials science (films, coatings)")

    return nothing
end

function analyze_phantom(phantom::AbstractArray, name::String,
                        grid_size::NTuple{3, Int}, resolution::NTuple{3, Float64})
    """Analyze and report phantom statistics with physical dimensions."""

    # Calculate material distribution
    background_voxels = count(real.(phantom) .≈ 1.0)
    object_voxels = count(real.(phantom) .> 1.1)
    total_voxels = prod(grid_size)

    volume_fraction = object_voxels / total_voxels

    # Physical volume calculations
    voxel_volume = prod(resolution)  # m³ per voxel
    object_volume = object_voxels * voxel_volume  # m³
    object_volume_μm3 = object_volume * 1e18  # μm³

    # Find permittivity range
    min_perm = minimum(real.(phantom))
    max_perm = maximum(real.(phantom))
    mean_perm = sum(real.(phantom)) / length(phantom)

    println("$name analysis:")
    println("  Grid size: $(grid_size)")
    println("  Object voxels: $object_voxels / $total_voxels ($(round(volume_fraction*100, digits=1))%)")
    println("  Object volume: $(round(object_volume_μm3, digits=3)) μm³")
    println("  Permittivity: $(round(min_perm, digits=3)) - $(round(max_perm, digits=3)) (mean: $(round(mean_perm, digits=3)))")

    # Refractive index information
    if max_perm > min_perm
        min_n = sqrt(min_perm)
        max_n = sqrt(max_perm)
        println("  Refractive index: $(round(min_n, digits=3)) - $(round(max_n, digits=3))")
    end

    # Contrast with respect to water
    water_perm = 1.33^2
    if max_perm > water_perm * 1.01  # If significantly different from water
        contrast = max_perm / water_perm - 1
        println("  Contrast vs water: $(round(contrast, digits=3))")
    end
end

function create_complex_phantom(grid_size::NTuple{3, Int})
    """Create a complex phantom with multiple objects representing a realistic scenario."""

    # Start with background (air/vacuum)
    phantom = ones(ComplexF64, grid_size)

    # Add glass substrate (bottom 1/4 of domain)
    z_substrate = div(grid_size[3], 4)
    phantom[:, :, 1:z_substrate] .= 1.52^2  # Borosilicate glass

    # Add water layer (middle portion)
    z_water_start = z_substrate + 1
    z_water_end = div(3 * grid_size[3], 4)
    phantom[:, :, z_water_start:z_water_end] .= 1.33^2  # Water

    # Add polystyrene beads in water layer
    # Create several beads at different positions
    centers = [
        (div(grid_size[1], 3), div(grid_size[2], 3), div(grid_size[3], 2)),
        (div(2*grid_size[1], 3), div(grid_size[2], 2), div(grid_size[3], 2)),
        (div(grid_size[1], 2), div(2*grid_size[2], 3), div(grid_size[3], 2))
    ]

    for (cx, cy, cz) in centers
        # Create bead region
        for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
            r = sqrt((i-cx)^2 + (j-cy)^2 + (k-cz)^2)
            if r ≤ 8.0  # 8-pixel radius
                phantom[i, j, k] = 1.59^2  # Polystyrene
            end
        end
    end

    return phantom
end

function create_cell_phantom(grid_size::NTuple{3, Int}, complexity::Symbol)
    """Create biological cell phantoms of varying complexity."""

    phantom = ones(ComplexF64, grid_size)  # Background (medium)

    # Cell boundary (large sphere)
    cell_center = (div(grid_size[1], 2), div(grid_size[2], 2), div(grid_size[3], 2))
    cell_radius = 20.0  # 20 pixels

    # Cytoplasm
    for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
        r = sqrt((i-cell_center[1])^2 + (j-cell_center[2])^2 + (k-cell_center[3])^2)
        if r ≤ cell_radius
            phantom[i, j, k] = 1.37^2  # Cytoplasm n=1.37
        end
    end

    # Nucleus (always present)
    nucleus_radius = 8.0  # 8 pixels
    for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
        r = sqrt((i-cell_center[1])^2 + (j-cell_center[2])^2 + (k-cell_center[3])^2)
        if r ≤ nucleus_radius
            phantom[i, j, k] = 1.39^2  # Nucleus n=1.39
        end
    end

    if complexity == :complex
        # Add organelles (mitochondria, etc.)
        organelle_centers = [
            (cell_center[1] + 10, cell_center[2], cell_center[3]),
            (cell_center[1] - 8, cell_center[2] + 6, cell_center[3]),
            (cell_center[1], cell_center[2] - 10, cell_center[3] + 4)
        ]

        for (ox, oy, oz) in organelle_centers
            for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
                r = sqrt((i-ox)^2 + (j-oy)^2 + (k-oz)^2)
                if r ≤ 3.0  # Small organelles
                    # Check if within cell boundary
                    r_cell = sqrt((i-cell_center[1])^2 + (j-cell_center[2])^2 + (k-cell_center[3])^2)
                    if r_cell ≤ cell_radius && r_cell > nucleus_radius
                        phantom[i, j, k] = 1.41^2  # Organelle n=1.41
                    end
                end
            end
        end
    end

    return phantom
end

# Run the example
phantom_gallery_example()
```

## Detailed Examples

### 1. Single Bead Phantoms

Create spherical particles for flow cytometry and particle sizing applications:

```julia
# Standard polystyrene calibration bead
bead_ps = phantom_bead(
    (128, 128, 64),      # Grid size
    [1.59^2],            # Polystyrene permittivity
    10.0                 # 10-pixel radius = 500 nm at 50 nm resolution
)

# Silica nanoparticle
bead_silica = phantom_bead(
    (128, 128, 64),
    [1.46^2],            # Silica permittivity
    6.0                  # 6-pixel radius = 300 nm
)

# Metallic nanoparticle (with loss)
n_gold = 0.47 + 2.4im    # Gold at 532 nm
bead_gold = phantom_bead(
    (128, 128, 64),
    [n_gold^2],          # Complex permittivity
    4.0                  # 4-pixel radius = 200 nm
)
```

### 2. Plate and Film Geometries

Model optical components and thin films:

```julia
# Cover glass for microscopy
cover_glass = phantom_plate(
    (128, 128, 64),
    [1.52^2],            # Borosilicate glass
    20.0,                # 20-pixel thickness = 1 μm
    :z                   # Horizontal orientation
)

# Anti-reflection coating
ar_coating = phantom_plate(
    (128, 128, 64),
    [1.22^2],            # MgF2 coating (n=1.22)
    3.0,                 # 3-pixel thickness = 150 nm
    :z                   # On top surface
)

# Vertical membrane
membrane = phantom_plate(
    (128, 128, 64),
    [1.35^2],            # Biological membrane
    2.0,                 # 2-pixel thickness = 100 nm
    :x                   # Vertical orientation
)
```

### 3. Cylindrical Structures

Model fibers, tubes, and biological filaments:

```julia
# Step-index optical fiber
fiber_core = phantom_cylinder(
    (128, 128, 64),
    [1.46^2],            # Core: silica
    8.0,                 # 8-pixel radius = 400 nm
    60.0,                # 60-pixel length = 6 μm
    :z                   # Along z-axis
)

# Cladding can be added separately or by modifying the phantom
# fiber_core[distance_from_axis > 8 && distance_from_axis < 12] = 1.44^2

# Microtubule
microtubule = phantom_cylinder(
    (128, 128, 64),
    [1.37^2],            # Protein density
    1.5,                 # 1.5-pixel radius = 75 nm (typical microtubule)
    40.0,                # 40-pixel length = 4 μm
    :z                   # Along z-axis
)
```

### 4. Biological Cell Models

Create realistic cell phantoms for biological imaging:

```julia
function create_biological_cell(grid_size, cell_type::Symbol)
    phantom = ones(ComplexF64, grid_size)
    center = (div(grid_size[1], 2), div(grid_size[2], 2), div(grid_size[3], 2))

    if cell_type == :red_blood_cell
        # Red blood cell (biconcave disc approximation)
        cell_radius = 15.0  # ~750 nm
        for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
            r_xy = sqrt((i-center[1])^2 + (j-center[2])^2)
            z_dist = abs(k - center[3])

            # Simplified biconcave shape
            if r_xy ≤ cell_radius && z_dist ≤ 3.0
                phantom[i, j, k] = 1.36^2  # RBC hemoglobin concentration
            end
        end

    elseif cell_type == :white_blood_cell
        # White blood cell with nucleus
        cell_radius = 18.0  # ~900 nm
        nucleus_radius = 8.0

        # Cell cytoplasm
        for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
            r = sqrt((i-center[1])^2 + (j-center[2])^2 + (k-center[3])^2)
            if r ≤ cell_radius
                phantom[i, j, k] = 1.37^2  # Cytoplasm
            end
            if r ≤ nucleus_radius
                phantom[i, j, k] = 1.39^2  # Nucleus
            end
        end
    end

    return phantom
end

# Usage
rbc = create_biological_cell((128, 128, 64), :red_blood_cell)
wbc = create_biological_cell((128, 128, 64), :white_blood_cell)
```

### 5. Multi-Object Phantoms

Combine multiple objects for realistic scenarios:

```julia
function create_flow_chamber(grid_size)
    """Simulate a microfluidic flow chamber with particles."""
    phantom = ones(ComplexF64, grid_size)

    # Glass walls (top and bottom)
    wall_thickness = 10
    phantom[:, :, 1:wall_thickness] .= 1.52^2         # Bottom wall
    phantom[:, :, end-wall_thickness+1:end] .= 1.52^2 # Top wall

    # Water in channel
    phantom[:, :, wall_thickness+1:end-wall_thickness] .= 1.33^2

    # Add flowing particles
    particle_positions = [
        (30, 40, 35), (60, 50, 35), (90, 45, 35),  # Row 1
        (40, 70, 35), (70, 75, 35), (100, 80, 35)  # Row 2
    ]

    for (px, py, pz) in particle_positions
        for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
            r = sqrt((i-px)^2 + (j-py)^2 + (k-pz)^2)
            if r ≤ 6.0  # 6-pixel radius particles
                phantom[i, j, k] = 1.59^2  # Polystyrene beads
            end
        end
    end

    return phantom
end

flow_chamber = create_flow_chamber((128, 128, 64))
```

## Material Properties Database

### Common Optical Materials

| Material | n | ε | Applications |
|----------|---|---|--------------|
| Air/Vacuum | 1.000 | 1.000 | Reference medium |
| Water | 1.333 | 1.777 | Biological imaging |
| Cytoplasm | 1.370 | 1.877 | Cell content |
| Cell nucleus | 1.390 | 1.932 | Dense cell regions |
| Silica (SiO₂) | 1.460 | 2.132 | Optical fibers, beads |
| Borosilicate | 1.520 | 2.310 | Cover glass, substrates |
| Polystyrene | 1.590 | 2.528 | Calibration beads |
| Titanium dioxide | 2.400 | 5.760 | High-index coatings |

### Material Selection Guidelines

**For biological imaging**:
- Background: Water (n = 1.333)
- Cells: Cytoplasm (n = 1.37), Nucleus (n = 1.39)
- Organelles: n = 1.41-1.43

**For flow cytometry**:
- Medium: Saline solution (n ≈ 1.335)
- Beads: Polystyrene (n = 1.59), Silica (n = 1.46)

**For materials characterization**:
- Substrates: Glass (n = 1.52), Silicon (n ≈ 4.0)
- Films: Various depending on application

## Phantom Validation and Analysis

### Statistical Analysis

```julia
function detailed_phantom_analysis(phantom, resolution)
    # Volume calculations
    total_voxels = length(phantom)
    object_voxels = count(real.(phantom) .> 1.1)
    volume_fraction = object_voxels / total_voxels

    # Physical dimensions
    voxel_volume = prod(resolution) * 1e18  # μm³
    object_volume = object_voxels * voxel_volume

    # Material contrast
    background_perm = mode(real.(phantom[real.(phantom) .< 1.5]))  # Most common low value
    max_contrast = maximum(real.(phantom)) / background_perm - 1

    # Spatial distribution
    center_of_mass = calculate_center_of_mass(phantom)

    println("Detailed phantom analysis:")
    println("  Total volume: $(round(total_voxels * voxel_volume, digits=2)) μm³")
    println("  Object volume: $(round(object_volume, digits=2)) μm³")
    println("  Volume fraction: $(round(volume_fraction * 100, digits=1))%")
    println("  Maximum contrast: $(round(max_contrast, digits=2))")
    println("  Center of mass: $(center_of_mass)")
end
```

### Visualization Recommendations

For phantom visualization, consider:

1. **Cross-sections**: Extract 2D slices for easy viewing
2. **Isosurfaces**: 3D rendering of material boundaries
3. **Histograms**: Permittivity distribution analysis
4. **Profiles**: Line plots through objects of interest

```julia
# Example visualization code
using Plots

# Extract central cross-section
z_center = div(size(phantom, 3), 2)
cross_section = real.(phantom[:, :, z_center])

# Plot with proper physical scaling
x_μm = (1:size(phantom, 1)) .* resolution[1] .* 1e6
y_μm = (1:size(phantom, 2)) .* resolution[2] .* 1e6

heatmap(x_μm, y_μm, cross_section',
        xlabel="x (μm)", ylabel="y (μm)",
        title="Phantom Cross-Section",
        aspect_ratio=:equal)
```

## Best Practices

### Phantom Design

1. **Resolution**: Use at least 10 points per wavelength in the medium
2. **Boundaries**: Leave adequate padding (>λ/2) from domain edges
3. **Materials**: Use realistic optical constants for your wavelength
4. **Contrast**: High contrast may require stricter convergence tolerances

### Performance Considerations

1. **Memory**: Large phantoms use significant RAM (16 bytes per voxel for ComplexF64)
2. **Precision**: Consider Float32 for memory-limited applications
3. **Geometry**: Simple shapes converge faster than complex structures

### Physical Realism

1. **Dispersion**: Material properties are wavelength-dependent
2. **Loss**: Include imaginary permittivity for absorbing materials
3. **Interfaces**: Smooth transitions often improve convergence
4. **Size**: Ensure objects are adequately sampled by the grid

This phantom gallery provides a comprehensive foundation for creating realistic material distributions in electromagnetic simulations using FrequencyMaxwell.jl.