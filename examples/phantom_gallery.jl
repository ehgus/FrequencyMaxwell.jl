"""
Phantom Gallery Example for FrequencyMaxwell.

This example demonstrates the various phantom generation capabilities
for creating test objects in electromagnetic simulations.
"""

using FrequencyMaxwell

function main()
    println("FrequencyMaxwell Phantom Gallery")
    println("=" ^ 35)

    # Common parameters
    grid_size = (128, 128, 64)
    resolution = (50e-9, 50e-9, 100e-9)  # 50nm x 50nm x 100nm voxels

    println("Grid parameters:")
    println("  Grid size: $grid_size")
    println("  Resolution: $(resolution .* 1e9) nm")
    println("  Domain size: $((grid_size .* resolution) .* 1e6) μm")
    println()

    # Example 1: Single polystyrene bead
    println("1. Single Polystyrene Bead")
    println("-" ^ 28)
    bead_single = phantom_bead(
        grid_size,
        [1.59^2],                  # Polystyrene n=1.59
        12.0,                      # 12-pixel radius (600 nm)
        num_bead = 1
    )

    analyze_phantom(bead_single, "Single bead", grid_size)
    println()

    # Example 2: Multiple beads with different materials
    println("2. Multiple Material Beads")
    println("-" ^ 27)
    bead_multi = phantom_bead(
        grid_size,
        [1.33^2, 1.59^2, 1.46^2], # Water, polystyrene, silica
        8.0,                       # 8-pixel radius (400 nm)
        num_bead = 5,              # 5 beads total
        bead_distance = 20.0       # 20-pixel spacing (1 μm)
    )

    analyze_phantom(bead_multi, "Multi-material beads", grid_size)
    println()

    # Example 3: Glass plate (cover slip)
    println("3. Glass Plate (Cover Slip)")
    println("-" ^ 28)
    plate_glass = phantom_plate(
        grid_size,
        [1.52^2],                  # Borosilicate glass n=1.52
        16.0,                      # 16-pixel thickness (1.6 μm)
        plate_normal = 3           # Z-normal (horizontal plate)
    )

    analyze_phantom(plate_glass, "Glass plate", grid_size)
    println()

    # Example 4: Thin film in Y direction
    println("4. Thin Film (Y-normal)")
    println("-" ^ 24)
    film_thin = phantom_plate(
        grid_size,
        [2.4^2],                   # High-index material n=2.4
        6.0,                       # 6-pixel thickness (300 nm)
        plate_normal = 2           # Y-normal (vertical film)
    )

    analyze_phantom(film_thin, "Thin film", grid_size)
    println()

    # Example 5: Cylindrical fiber
    println("5. Cylindrical Optical Fiber")
    println("-" ^ 31)
    fiber = phantom_cylinder(
        grid_size,
        1.46^2,                    # Silica core n=1.46
        8.0,                       # 8-pixel radius (400 nm)
        40.0,                      # 40-pixel height (4 μm)
        axis = 3                   # Z-axis (vertical fiber)
    )

    analyze_phantom(fiber, "Optical fiber", grid_size)
    println()

    # Example 6: Comparison of materials
    println("6. Material Database Examples")
    println("-" ^ 30)

    materials = [
        ("Air/Vacuum", 1.0^2),
        ("Water", 1.33^2),
        ("Cytoplasm", 1.37^2),
        ("Silica", 1.46^2),
        ("Borosilicate", 1.52^2),
        ("Polystyrene", 1.59^2),
        ("Titanium Dioxide", 2.4^2)
    ]

    println("Common optical materials:")
    for (name, perm) in materials
        n = sqrt(real(perm))
        contrast = perm / 1.33^2 - 1  # Contrast vs water
        println("  $name: n = $(round(n, digits=3)), ε = $(round(real(perm), digits=3)), Δε = $(round(real(contrast), digits=3))")
    end
    println()

    # Example 7: Complex phantom with multiple objects
    println("7. Complex Multi-Object Phantom")
    println("-" ^ 33)

    # Start with background
    complex_phantom = ones(ComplexF64, grid_size)

    # Add glass substrate
    substrate = phantom_plate(grid_size, [1.52^2], 20.0, plate_normal = 3)
    complex_phantom[real.(substrate) .> 1.1] = 1.52^2

    # Add multiple beads on substrate
    for i in 1:3
        bead_pos = phantom_bead(
            grid_size,
            [1.59^2],
            6.0,
            num_bead = 1
        )
        # Shift bead position (this is a simplified approach)
        complex_phantom[real.(bead_pos) .> 1.1] = 1.59^2
    end

    analyze_phantom(complex_phantom, "Complex multi-object", grid_size)
    println()

    println("Phantom Gallery Summary:")
    println("- All phantoms use realistic optical material properties")
    println("- Phantoms are suitable for electromagnetic scattering simulations")
    println("- Different geometries demonstrate various applications:")
    println("  * Beads: Cell imaging, flow cytometry")
    println("  * Plates: Cover slips, optical windows")
    println("  * Cylinders: Optical fibers, biological filaments")
    println("  * Complex: Realistic experimental scenarios")

    return nothing
end

function analyze_phantom(phantom::AbstractArray, name::String, grid_size::NTuple{3, Int})
    """Analyze and report phantom statistics."""

    # Calculate material distribution
    background_voxels = count(real.(phantom) .≈ 1.0)
    object_voxels = count(real.(phantom) .> 1.1)
    total_voxels = prod(grid_size)

    volume_fraction = object_voxels / total_voxels

    # Find permittivity range
    min_perm = minimum(real.(phantom))
    max_perm = maximum(real.(phantom))
    mean_perm = sum(real.(phantom)) / length(phantom)

    println("$name statistics:")
    println("  Size: $(size(phantom))")
    println("  Background voxels: $background_voxels ($(round(background_voxels/total_voxels*100, digits=1))%)")
    println("  Object voxels: $object_voxels ($(round(volume_fraction*100, digits=1))%)")
    println("  Permittivity range: $(round(min_perm, digits=3)) - $(round(max_perm, digits=3))")
    println("  Mean permittivity: $(round(mean_perm, digits=3))")

    # Calculate refractive index info if meaningful
    if max_perm > min_perm
        min_n = sqrt(min_perm)
        max_n = sqrt(max_perm)
        println("  Refractive index range: $(round(min_n, digits=3)) - $(round(max_n, digits=3))")
    end
end

# Run the example if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
