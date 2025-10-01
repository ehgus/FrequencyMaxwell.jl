"""
Two-Beam Interference Example for Multi-Source Electromagnetic Solving

This example demonstrates multi-source electromagnetic simulation capability by reproducing
the MATLAB two-beam examples (forward_two_beam_plain_medium.m and forward_two_beam_SiO2_5um_in_water.m).

The simulation shows coherent interference between two plane wave sources with opposite
horizontal k-vectors, creating characteristic interference fringes in the simulation domain.

# Key Features Demonstrated:
- Multi-source coherent electromagnetic simulation
- Interference pattern generation from coherent superposition
- FrequencyMaxwell's type-safe multi-source API
- Validation against MATLAB reference implementation

# Physical Setup:
- Two plane waves with ±3rd order horizontal angles
- 532 nm wavelength in water background (n=1.333)
- X-polarized electromagnetic fields
- Coherent superposition producing interference fringes
"""

using LinearAlgebra
using FrequencyMaxwell
using LinearSolve
using WGLMakie

"""
    test_two_beam_interference_homogeneous()

Test two-beam interference in homogeneous medium (water).

This reproduces the MATLAB forward_two_beam_plain_medium.m example.
Expected result: sinusoidal interference pattern with high contrast fringes.
"""
function test_two_beam_interference_homogeneous()
    println("=== Two-Beam Interference Test (Homogeneous Medium) ===")

    # Configure solver (matching MATLAB parameters)
    # Define boundary conditions: periodic in X,Y and absorbing in Z
    bc_absorbing_z = AbsorbingBoundaryCondition(
        thickness = 300e-9,              # 300nm padding (6 pixels * 50nm)
        attenuation_thickness = 300e-9,  # 300nm attenuation layer
        sharpness = 1.0,                 # Sharp attenuation
        profile = TanhProfile            # Smooth tanh profile
    )

    solver = ConvergentBornSolver(
        wavelength = 532e-9,           # 532 nm
        permittivity_bg = 1.333^2,     # Water background (n=1.333)
        resolution = (50e-9, 50e-9, 50e-9),  # 50 nm isotropic resolution
        grid_size = (201, 201, 191),   # Matching MATLAB grid
        boundary_conditions = (        # Periodic in XY, absorbing in Z
            PeriodicBoundaryCondition(),
            PeriodicBoundaryCondition(),
            bc_absorbing_z
        ),
        iterations_max = -1,           # Auto-determine iterations
        tolerance = 1e-6,
        linear_solver = KrylovJL_BICGSTAB()
    )

    println("✓ Solver configured with:")
    println("  Wavelength: $(solver.wavelength*1e9) nm")
    println("  Background permittivity: $(solver.permittivity_bg)")
    println("  Grid size: $(solver.grid_size)")
    println("  Domain size: $(domain_size(solver) .* 1e6) μm")

    # Create two-beam sources (matching MATLAB illum_order = 3)
    illum_order = 3
    ky = 2π * illum_order / (solver.grid_size[2] * solver.resolution[2])

    # Calculate propagation angle
    k_bg = wavenumber_background(solver)
    angle = asin(ky / k_bg)

    println("✓ Two-beam configuration:")
    println("  Illumination order: $(illum_order)")
    println("  Horizontal k-vector: ±$(ky) rad/m")
    println("  Beam angle: ±$(rad2deg(angle))°")

    # Create plane wave sources with opposite horizontal angles
    source1 = PlaneWaveSource(
        wavelength = solver.wavelength,
        polarization = [1.0, 0.0, 0.0],  # X-polarized
        k_vector = [0.0, sin(angle), cos(angle)],  # +angle beam
        amplitude = 1.0
    )

    source2 = PlaneWaveSource(
        wavelength = solver.wavelength,
        polarization = [1.0, 0.0, 0.0],  # X-polarized
        k_vector = [0.0, -sin(angle), cos(angle)], # -angle beam
        amplitude = 1.0
    )

    sources = [source1, source2]

    # Validate source configuration
    println("✓ Sources created:")
    for (i, src) in enumerate(sources)
        println("  Source $(i): k = $(src.k_vector), λ = $(src.wavelength*1e9) nm")
    end

    # Homogeneous medium (water everywhere)
    permittivity = fill(solver.permittivity_bg, solver.grid_size)

    println("✓ Homogeneous medium: ε = $(solver.permittivity_bg)")

    # Solve multi-source electromagnetic problem
    println("\n--- Solving Multi-Source CBS Problem ---")
    t_start = time()
    EMfield = solve(solver, sources, permittivity)
    t_solve = time() - t_start

    println("✓ Multi-source solve completed in $(t_solve) seconds")
    println("  E-field size: $(size(EMfield.E))")
    println("  H-field size: $(size(EMfield.H))")

    # Analyze interference pattern
    intensity = sum(abs2.(EMfield.E), dims = 4)[:, :, :, 1]

    println("\n--- Analyzing Interference Pattern ---")

    # Extract central Y cross-section for fringe analysis
    center_x = div(size(intensity, 1), 2) + 1
    center_z = div(size(intensity, 3), 2) + 1
    central_y_profile = intensity[center_x, :, center_z]

    # Calculate fringe statistics
    max_intensity = maximum(central_y_profile)
    min_intensity = minimum(central_y_profile)
    contrast = (max_intensity - min_intensity) / (max_intensity + min_intensity)

    println("✓ Interference pattern analysis:")
    println("  Max intensity: $(max_intensity)")
    println("  Min intensity: $(min_intensity)")
    println("  Fringe contrast: $(contrast)")
    println("  Expected contrast: ~1.0 for perfect interference")

    # Count interference fringes
    # For illum_order = 3, expect ~6 fringes across Y dimension
    y_coords = (0:(solver.grid_size[2] - 1)) .* solver.resolution[2]
    expected_fringes = 2 * illum_order

    # Simple fringe counting by detecting maxima
    fringe_count = count_intensity_maxima(central_y_profile)

    println("✓ Fringe analysis:")
    println("  Detected fringes: $(fringe_count)")
    println("  Expected fringes: $(expected_fringes)")
    println("  Fringe spacing: $(solver.grid_size[2] * solver.resolution[2] / expected_fringes * 1e6) μm")

    # Validation checks
    success = true

    if contrast < 0.8
        println("❌ WARNING: Low contrast $(contrast), expected >0.8")
        success = false
    else
        println("✓ High contrast interference achieved")
    end

    if abs(fringe_count - expected_fringes) > 2
        println("❌ WARNING: Fringe count $(fringe_count) differs from expected $(expected_fringes)")
        success = false
    else
        println("✓ Fringe count matches expected pattern")
    end

    # Check field magnitudes are reasonable
    E_max = maximum(abs.(EMfield.E))
    H_max = maximum(abs.(EMfield.H))

    if E_max < 1e-3 || E_max > 1e3
        println("❌ WARNING: E-field magnitude $(E_max) seems unrealistic")
        success = false
    else
        println("✓ E-field magnitude reasonable: $(E_max)")
    end

    println("\n=== Multi-Source Interference Test Results ===")
    if success
        println("🎉 SUCCESS: Multi-source interference working correctly!")
        println("   - Coherent superposition implemented properly")
        println("   - High-contrast interference fringes generated")
        println("   - Fringe pattern matches theoretical expectation")
        println("   - Field magnitudes are physically reasonable")
    else
        println("⚠️  PARTIAL: Some validation checks failed, review above warnings")
    end

    return EMfield, intensity, success
end

"""
    count_intensity_maxima(profile)

Count the number of intensity maxima in a 1D intensity profile.
Used for fringe counting in interference patterns.
"""
function count_intensity_maxima(profile::AbstractVector)
    maxima_count = 0
    n = length(profile)

    for i in 2:(n - 1)
        if profile[i] > profile[i - 1] && profile[i] > profile[i + 1]
            # Additional check: make sure it's a significant maximum
            if profile[i] > 0.1 * maximum(profile)
                maxima_count += 1
            end
        end
    end

    return maxima_count
end

"""
    create_two_beam_sources(wavelength, angle_deg; polarization=[1,0,0])

Convenience function to create two plane wave sources with opposite horizontal angles.
Matches the MATLAB PlaneSource configuration pattern.
"""
function create_two_beam_sources(wavelength, angle_deg; polarization = [1.0, 0.0, 0.0])
    angle_rad = deg2rad(angle_deg)

    source1 = PlaneWaveSource(
        wavelength = wavelength,
        polarization = polarization,
        k_vector = [0.0, sin(angle_rad), cos(angle_rad)]
    )

    source2 = PlaneWaveSource(
        wavelength = wavelength,
        polarization = polarization,
        k_vector = [0.0, -sin(angle_rad), cos(angle_rad)]
    )

    return [source1, source2]
end

# Run the test if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    try
        EMfield, intensity, success = test_two_beam_interference_homogeneous()
        heatmap(angle.(EMfield.E[div(size(EMfield.E, 1), 2), :, :, 1]))
        heatmap(intensity[div(size(intensity, 1), 2), :, :])
        exit(success ? 0 : 1)
    catch e
        println("❌ MULTI-SOURCE TEST FAILED: $e")
        println("\nStack trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        exit(1)
    end
end
