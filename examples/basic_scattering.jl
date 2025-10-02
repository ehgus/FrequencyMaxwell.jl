"""
Basic electromagnetic scattering example using FrequencyMaxwell.

This example demonstrates the basic usage of FrequencyMaxwell for simulating
electromagnetic scattering from a simple spherical bead phantom.
"""

using FrequencyMaxwell
using LinearSolve
using WGLMakie

function main()
    println("FrequencyMaxwell Basic Scattering Example")
    println("=" ^ 45)

    # Configure and create the electromagnetic solver with optimized CBS parameters
    println("Creating solver with enhanced API...")

    # Define boundary conditions: periodic in X,Y and absorbing in Z
    bc_absorbing_z = AbsorbingBoundaryCondition(
        thickness = 3.0e-6,              # 3 μm padding in Z
        attenuation_thickness = 3.0e-6,  # 3 μm attenuation layer
        sharpness = 1.0,                 # Sharp attenuation
        profile = TanhProfile            # Smooth tanh profile (default)
    )

    solver = ConvergentBornSolver(
        permittivity_bg = 1.333^2, # Water background (n=1.333)
        resolution = (50e-9, 50e-9, 50e-9),  # 50 nm isotropic resolution
        grid_size = (128, 128, 64),           # 128×128×64 grid
        boundary_conditions = (               # Periodic in XY, absorbing in Z
            PeriodicBoundaryCondition(),
            PeriodicBoundaryCondition(),
            bc_absorbing_z
        ),
        iterations_max = -1,       # Auto-calculate optimal iterations
        tolerance = 1e-2,          # Convergence tolerance
        linear_solver = KrylovJL_GMRES()  # LinearSolver object
    )

    println("Configuration:")
    println("  Background n: $(sqrt(solver.permittivity_bg))")
    println("  Grid size: $(solver.grid_size)")
    println("  Resolution: $(solver.resolution .* 1e9) nm")
    println("  Domain size: $(domain_size(solver) .* 1e6) μm")

    # Define the incident plane wave source
    println("Setting up plane wave source...")
    source = PlaneWaveSource(
        wavelength = 532e-9,      # 532 nm (green laser)
        polarization = [1.0, 0.0, 0.0],    # X-polarized
        k_vector = [0.0, 0.0, 1.0],        # Propagating in +Z
        amplitude = 1.0                     # 1 V/m amplitude
    )

    println("Source configuration:")
    println("  Wavelength: $(source_wavelength(source) * 1e9) nm")
    println("  Polarization: $(source.polarization)")
    println("  Propagation: $(source.k_vector)")

    # Visualize the plane wave source
    incident_field = generate_incident_field(source, solver)
    heatmap(angle.(incident_field.E[div(size(incident_field.E, 1), 2), :, :, 1]))

    # Create a spherical bead phantom (smaller for more realistic scattering)
    println("\nGenerating phantom...")
    bead_radius_pixels = 10.0  # 500 nm radius (10 pixels * 50 nm)
    phantom = phantom_bead(
        solver.grid_size,
        [1.46^2],                  # SiO2 bead (n=1.46, realistic optical material)
        bead_radius_pixels        # 500 nm radius
    )

    # Calculate phantom statistics
    bead_volume = count(real.(phantom) .> 1.1)  # Voxels with elevated permittivity
    total_volume = prod(solver.grid_size)
    volume_fraction = bead_volume / total_volume

    println("Phantom statistics:")
    println("  Bead material: n = $(sqrt(1.46^2)) (silica)")
    println("  Bead radius: $(bead_radius_pixels * solver.resolution[1] * 1e9) nm")
    println("  Volume fraction: $(round(volume_fraction * 100, digits=2))%")

    heatmap(abs.(phantom[div(size(phantom, 1), 2), :, :]))

    # Solve the electromagnetic scattering problem
    println("\nSolving electromagnetic scattering...")
    println("Note: This is a placeholder implementation")

    EMfield = solve(solver, source, phantom)

    println("Solution completed successfully!")
    println("Field properties:")
    println("  Electric field size: $(size(EMfield.E))")
    println("  Magnetic field size: $(size(EMfield.H))")
    println("  Total field energy: $(field_energy(EMfield)) J")

    # Calculate field intensity
    intensity = field_intensity(EMfield)
    max_intensity = maximum(intensity)
    enhancement = max_intensity / 1.0  # Relative to incident intensity

    heatmap(intensity[div(size(intensity, 1), 2), :, :])
    heatmap(angle.(EMfield.E[div(size(EMfield.E, 1), 2), :, :, 1]))
    println("  Maximum intensity: $(max_intensity) (V/m)²")
    println("  Enhancement factor: $(round(enhancement, digits=2))")

    # Extract central plane for analysis
    central_plane = extract_plane(EMfield, 3, div(solver.grid_size[3], 2))
    plane_intensity = field_intensity(central_plane)

    println("Central plane analysis:")
    println("  Plane size: $(size(plane_intensity))")
    println("  Mean intensity: $(round(sum(plane_intensity)/length(plane_intensity), digits=4))")

    println("\nExample completed successfully!")

    return nothing
end

# Run the example if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
