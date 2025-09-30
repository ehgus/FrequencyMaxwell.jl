"""
Forward solver integration tests - migrated from forward_SiO2_5um_in_water.jl example
Tests complete forward solver workflows including realistic scattering problems.
"""

@testset "Forward Solver Integration" begin
    @testset "SiO2 Bead in Water - Complete Workflow" begin
        # Replicate the original example with smaller grid for testing
        radius_um = 0.5  # 0.5 micron radius
        solver = ConvergentBornSolver(
            wavelength = 532e-9,
            permittivity_bg = 1.333^2,  # Water
            resolution = (50e-9, 50e-9, 50e-9),  # 50 nm resolution
            grid_size = (101, 101, 96),  # Smaller than original 201x201x191
            boundary_thickness = (0, 0, 6),
            field_attenuation = (0, 0, 6)
        )

        # Create SiO2 bead phantom
        permittivity_SiO2 = 1.4607^2
        radius_pixels = round(Int, radius_um * 1e-6 / solver.resolution[3])
        phantom = create_spherical_phantom(
            solver.grid_size,
            solver.permittivity_bg,
            permittivity_SiO2,
            radius_pixels
        )

        # Create plane wave source
        source = PlaneWaveSource(
            solver.wavelength,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )

        # Set permittivity and solve forward problem
        set_permittivity!(solver, phantom)

        # Test solver state after permittivity setting
        @test !isnothing(solver.permittivity)
        @test size(solver.permittivity) == solver.grid_size

        # Solve with limited Born iterations for testing
        set_Born_max!(solver, 10)  # Reduced from default for speed
        E_field, H_field = solve(solver, [source])

        @test size(E_field) == (solver.grid_size..., 3, 1)  # Last dimension for sources
        @test size(H_field) == (solver.grid_size..., 3, 1)

        # Test field properties
        @test !any(isnan.(E_field))
        @test !any(isnan.(H_field))
        @test !any(isinf.(E_field))
        @test !any(isinf.(H_field))

        # Calculate intensity
        intensity = sum(abs2.(E_field[:, :, :, :, 1]), dims = 4)

        # Test intensity properties
        @test all(intensity .≥ 0)
        @test maximum(intensity) > minimum(intensity)  # Should have contrast

        # Test that scattering enhances/reduces intensity in different regions
        background_intensity = intensity[1, 1, 1]  # Far from bead
        center_intensity = intensity[div(solver.grid_size[1], 2) + 1, div(
            solver.grid_size[2], 2) + 1, div(solver.grid_size[3], 2) + 1]

        @test background_intensity > 0
        @test center_intensity != background_intensity  # Scattering effect
    end

    @testset "Convergence Properties" begin
        # Test Born series convergence with small perturbation
        solver = ConvergentBornSolver(
            wavelength = 500e-9,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32)
        )

        # Small perturbation for good convergence
        phantom = ones(ComplexF64, solver.grid_size) * solver.permittivity_bg
        phantom[16:18, 16:18, 16:18] .= 1.01  # 1% permittivity increase

        source = PlaneWaveSource(
            solver.wavelength,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )
        set_permittivity!(solver, phantom)

        # Test different Born iteration counts
        results = []
        for Born_max in [1, 3, 5, 10]
            set_Born_max!(solver, Born_max)
            E_field, _ = solve(solver, [source])
            intensity = sum(abs2.(E_field[:, :, :, :, 1]), dims = 4)
            push!(results, intensity)
        end

        # Test that results converge (later iterations closer to each other)
        diff_1_3 = maximum(abs.(results[1] - results[2]))
        diff_5_10 = maximum(abs.(results[3] - results[4]))
        @test diff_5_10 < diff_1_3  # Better convergence with more iterations
    end

    @testset "Material Contrast Effects" begin
        solver = ConvergentBornSolver(
            wavelength = 633e-9,
            permittivity_bg = 1.33^2,  # Water
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32)
        )

        source = PlaneWaveSource(
            solver.wavelength,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )
        set_Born_max!(solver, 5)

        # Test different material contrasts
        contrasts = [1.01, 1.1, 1.5, 2.0]  # Increasing contrast
        scattering_strengths = []

        for contrast in contrasts
            phantom = ones(ComplexF64, solver.grid_size) * solver.permittivity_bg
            phantom[14:18, 14:18, 14:18] .= solver.permittivity_bg * contrast

            set_permittivity!(solver, phantom)
            E_field, _ = solve(solver, [source])

            intensity = sum(abs2.(E_field[:, :, :, :, 1]), dims = 4)
            background = mean(intensity[1:5, 1:5, 1:5])
            scattered = mean(intensity[14:18, 14:18, 14:18])

            scattering_strength = abs(scattered - background) / background
            push!(scattering_strengths, scattering_strength)
        end

        # Test that higher contrast produces stronger scattering
        @test all(diff(scattering_strengths) .> 0)  # Monotonically increasing
    end

    @testset "Lossy Materials" begin
        solver = ConvergentBornSolver(
            wavelength = 500e-9,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32)
        )

        # Test with lossy material (complex permittivity)
        phantom = ones(ComplexF64, solver.grid_size)
        phantom[12:20, 12:20, 12:20] .= 2.0 + 0.1im  # Lossy material

        source = PlaneWaveSource(
            solver.wavelength,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )
        set_permittivity!(solver, phantom)
        set_Born_max!(solver, 5)

        E_field, H_field = solve(solver, [source])

        # Test that fields are computed for lossy materials
        @test !any(isnan.(E_field))
        @test !any(isnan.(H_field))

        # Test absorption effects
        intensity = sum(abs2.(E_field[:, :, :, :, 1]), dims = 4)

        # Compare with lossless case
        phantom_lossless = ones(ComplexF64, solver.grid_size)
        phantom_lossless[12:20, 12:20, 12:20] .= 2.0  # Same real part

        set_permittivity!(solver, phantom_lossless)
        E_lossless, _ = solve(solver, [source])
        intensity_lossless = sum(abs2.(E_lossless[:, :, :, :, 1]), dims = 4)

        # Lossy material should generally reduce intensity
        mean_intensity_lossy = mean(intensity)
        mean_intensity_lossless = mean(intensity_lossless)
        @test mean_intensity_lossy <= mean_intensity_lossless
    end

    @testset "Field Boundary Conditions" begin
        solver = ConvergentBornSolver(
            wavelength = 500e-9,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32),
            boundary_thickness = (2, 2, 2),
            field_attenuation = (2, 2, 2)
        )

        # Simple scatterer
        phantom = ones(ComplexF64, solver.grid_size)
        phantom[15:17, 15:17, 15:17] .= 1.5

        source = PlaneWaveSource(
            solver.wavelength,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )
        set_permittivity!(solver, phantom)
        set_Born_max!(solver, 3)

        E_field, _ = solve(solver, [source])

        # Test that boundary attenuation is applied
        edge_intensity = mean(abs2.(E_field[1:2, :, :, 1, 1]))
        center_intensity = mean(abs2.(E_field[15:17, 15:17, 15:17, 1, 1]))

        # Edge should have lower intensity due to attenuation
        @test edge_intensity < center_intensity
    end
end
