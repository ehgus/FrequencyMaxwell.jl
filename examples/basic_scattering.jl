"""
Basic electromagnetic scattering example using FrequencyMaxwell.

This example demonstrates the basic usage of FrequencyMaxwell for simulating
electromagnetic scattering from a simple spherical bead phantom.
"""

using FrequencyMaxwell
using WGLMakie

function main()
    println("FrequencyMaxwell Basic Scattering Example")
    println("=" ^ 45)
    
    # Configure the electromagnetic solver
    println("Setting up solver configuration...")
    config = ConvergentBornConfig(
        wavelength = 500e-9,      # 500 nm (green light)
        permittivity_bg = 1.33^2, # Water background (n=1.33)
        resolution = (50e-9, 50e-9, 50e-9),  # 50 nm isotropic resolution
        grid_size = (128, 128, 64),          # 128×128×64 grid
        tolerance = 1e-6          # Convergence tolerance
    )
    
    println("Configuration:")
    println("  Wavelength: $(config.wavelength * 1e9) nm")
    println("  Background n: $(sqrt(config.permittivity_bg))")
    println("  Grid size: $(config.grid_size)")
    println("  Resolution: $(config.resolution .* 1e9) nm")
    println("  Domain size: $(domain_size(config) .* 1e6) μm")
    
    # Create the electromagnetic solver
    println("\nCreating solver...")
    solver = ConvergentBornSolver(config)
    
    # Define the incident plane wave source
    println("Setting up plane wave source...")
    source = PlaneWaveSource(
        wavelength = config.wavelength,
        polarization = [1.0, 0.0, 0.0],    # X-polarized
        k_vector = [0.0, 0.0, 1.0],        # Propagating in +Z
        amplitude = 1.0                     # 1 V/m amplitude
    )
    
    println("Source configuration:")
    println("  Wavelength: $(source_wavelength(source) * 1e9) nm")
    println("  Polarization: $(source.polarization)")
    println("  Propagation: $(source.k_vector)")
    println("  Power density: $(source_power(source)) W/m²")
    
    # Create a spherical bead phantom
    println("\nGenerating phantom...")
    phantom = phantom_bead(
        config.grid_size,
        [1.5^2],                   # Polystyrene bead (n=1.5)
        16.0,                      # 16-pixel radius (800 nm)
    )
    
    # Calculate phantom statistics
    bead_volume = count(real.(phantom) .> 1.1)  # Voxels with elevated permittivity
    total_volume = prod(config.grid_size)
    volume_fraction = bead_volume / total_volume
    
    println("Phantom statistics:")
    println("  Bead material: n = $(sqrt(1.5^2)) (polystyrene)")
    println("  Bead radius: $(16 * config.resolution[1] * 1e9) nm")
    println("  Volume fraction: $(round(volume_fraction * 100, digits=2))%")

    heatmap(abs.(phantom[:,:,div(config.grid_size[3],2)]))
    
    # Solve the electromagnetic scattering problem
    println("\nSolving electromagnetic scattering...")
    println("Note: This is a placeholder implementation")
    
    E_field, H_field = solve(solver, source, phantom)
    
    # Wrap fields in structured container
    fields = ElectromagneticField(E_field, H_field, config.grid_size, 
                                config.resolution, config.wavelength)
    
    println("Solution completed successfully!")
    println("Field properties:")
    println("  Electric field size: $(size(fields.E))")
    println("  Magnetic field size: $(size(fields.H))")
    println("  Total field energy: $(field_energy(fields)) J")
    
    # Calculate field intensity
    intensity = field_intensity(fields)
    max_intensity = maximum(intensity)
    enhancement = max_intensity / 1.0  # Relative to incident intensity
    
    heatmap(intensity[:,:,div(config.grid_size[3],2)])
    heatmap(angle.(fields.E[:,:,div(config.grid_size[3],2),1]))
    println("  Maximum intensity: $(max_intensity) (V/m)²")
    println("  Enhancement factor: $(round(enhancement, digits=2))")
    
    # Extract central plane for analysis
    central_plane = extract_plane(fields, 3, div(config.grid_size[3], 2))
    plane_intensity = field_intensity(central_plane)
    
    println("Central plane analysis:")
    println("  Plane size: $(size(plane_intensity))")
    println("  Mean intensity: $(round(sum(plane_intensity)/length(plane_intensity), digits=4))")
    
    println("\nExample completed successfully!")
    println("This demonstrates the basic FrequencyMaxwell workflow:")
    println("1. Configure solver parameters")
    println("2. Create electromagnetic solver")
    println("3. Define incident field source")
    println("4. Generate object phantom")
    println("5. Solve scattering problem")
    println("6. Analyze results")
    
    return nothing
end

# Run the example if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
