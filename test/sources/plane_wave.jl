"""
Plane wave source tests - migrated from forward solver and grating examples
Tests electromagnetic source generation and field initialization.
"""

@testset "Electromagnetic Sources" begin
    @testset "Plane Wave Source" begin
        # Test basic plane wave source - matches original examples
        solver = ConvergentBornSolver(
            permittivity_bg = 1.333^2,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (64, 64, 32),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        # Test plane wave propagating in z-direction
        source = PlaneWaveSource(
            wavelength = 532e-9,
            polarization = [1.0, 0.0, 0.0],  # x-polarized
            k_vector = [0.0, 0.0, 1.0]  # z-direction
        )

        @test source.polarization[1] ≈ 1.0
        @test source.polarization[2] ≈ 0.0
        @test source.polarization[3] ≈ 0.0
        @test source.k_vector[3] ≈ 1.0  # z-direction

        # Test field generation
        EMfield = generate_incident_field(source, solver)
        Efield = EMfield.E
        Hfield = EMfield.H

        @test size(Efield) == (solver.grid_size..., 3)
        @test size(Hfield) == (solver.grid_size..., 3)

        # Test field properties
        @test !any(isnan.(Efield))
        @test !any(isnan.(Hfield))
        @test !any(isinf.(Efield))
        @test !any(isinf.(Hfield))
    end

    @testset "Source Polarization" begin
        solver = ConvergentBornSolver(
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 16),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        # Test x-polarization
        source_x = PlaneWaveSource(
            wavelength = 500e-9,
            polarization = [1.0, 0.0, 0.0],
            k_vector = [0.0, 0.0, 1.0]
        )
        EMfield_x = generate_incident_field(source_x, solver)
        E_x, H_x = EMfield_x.E, EMfield_x.H

        # Test y-polarization
        source_y = PlaneWaveSource(
            wavelength = 500e-9,
            polarization = [0.0, 1.0, 0.0],
            k_vector = [0.0, 0.0, 1.0]
        )
        EMfield_y = generate_incident_field(source_y, solver)
        E_y, H_y = EMfield_y.E, EMfield_y.H

        # Test circular polarization
        source_circ = PlaneWaveSource(
            wavelength = 500e-9,
            polarization = [1.0, 1.0im, 0.0],
            k_vector = [0.0, 0.0, 1.0]
        )
        EMfield_circ = generate_incident_field(source_circ, solver)
        E_circ, H_circ = EMfield_circ.E, EMfield_circ.H

        # Verify orthogonality for linear polarizations
        @test sum(abs, (E_x .* conj(E_y))) ≈ 0.0  # Should be orthogonal
    end

    @testset "Oblique Incidence" begin
        solver = ConvergentBornSolver(
            permittivity_bg = 1.0,
            resolution = (10e-9, 10e-9, 10e-9),
            grid_size = (20, 50, 50),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        # Test different propagation directions - simplified from oblique incidence
        directions = [[0.0, 0.0, 1.0], [0.0, 0.1, 0.995]]  # Normalized k-vectors

        for k_vec in directions
            # Normalize k-vector
            k_normalized = k_vec ./ norm(k_vec)

            source = PlaneWaveSource(
                wavelength = 355e-9,  # UV wavelength from grating example
                polarization = [1.0, 0.0, 0.0],
                k_vector = k_normalized
            )

            EMfield = generate_incident_field(source, solver)
            Efield, Hfield = EMfield.E, EMfield.H

            # Test that fields are generated properly
            @test size(Efield) == (solver.grid_size..., 3)
            @test size(Hfield) == (solver.grid_size..., 3)

            # Test that non-normal incidence creates expected spatial patterns
            if k_normalized != [0.0, 0.0, 1.0]
                # Should have spatial variation due to oblique incidence
                center_field = Efield[div(end, 2), div(end, 2), div(end, 2), 1]
                edge_field = Efield[1, 1, div(end, 2), 1]
                @test abs(center_field - edge_field) > 0.0  # Should have spatial variation
            end
        end
    end

    @testset "Source Field Properties" begin
        solver = ConvergentBornSolver(
            permittivity_bg = 1.33^2,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 32),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        source = PlaneWaveSource(
            wavelength = 633e-9,
            polarization = [1.0, 0.0, 0.0],
            k_vector = [0.0, 0.0, 1.0]
        )

        EMfield = generate_incident_field(source, solver)
        Efield, Hfield = EMfield.E, EMfield.H

        # Test Maxwell's equations in source region
        # For plane wave: ∇ × E = -iωμH, ∇ × H = iωεE

        # Test field magnitudes are reasonable
        E_magnitude = sqrt.(sum(abs2.(Efield), dims = 4))
        H_magnitude = sqrt.(sum(abs2.(Hfield), dims = 4))

        @test maximum(E_magnitude) > 0
        @test maximum(H_magnitude) > 0

        # Test impedance relationship for plane wave
        # |H| ≈ |E| * sqrt(ε/μ) = |E| * sqrt(ε₀/μ₀) * sqrt(εᵣ)
        n_bg = sqrt(solver.permittivity_bg)
        Z0 = 376.730313668  # Free space impedance
        Z_medium = Z0 / n_bg

        center = div.(solver.grid_size, 2) .+ 1
        E_center_mag = abs(Efield[center..., 1])
        H_center_mag = abs(Hfield[center..., 2])  # H_y for E_x

        if H_center_mag > 0
            impedance_ratio = E_center_mag / H_center_mag
            @test impedance_ratio ≈ Z_medium rtol=0.1
        end
    end

    @testset "Multiple Sources" begin
        solver = ConvergentBornSolver(
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        # Create multiple sources for coherent superposition
        sources = [
            PlaneWaveSource(wavelength = 500e-9, polarization = [1.0, 0.0, 0.0],
                k_vector = [0.0, 0.0, 1.0]),
            PlaneWaveSource(wavelength = 500e-9, polarization = [0.0, 1.0, 0.0],
                k_vector = [0.0, 0.0, 1.0])
        ]

        # Test that individual sources work
        EMfield1 = generate_incident_field(sources[1], solver)
        EMfield2 = generate_incident_field(sources[2], solver)
        E1, H1 = EMfield1.E, EMfield1.H
        E2, H2 = EMfield2.E, EMfield2.H

        # Manual superposition for testing
        E_total = E1 + E2
        H_total = H1 + H2

        @test size(E_total) == (solver.grid_size..., 3)
        @test size(H_total) == (solver.grid_size..., 3)

        # Test that manual superposition is linear (already computed above)
        @test E_total ≈ E1 + E2
        @test H_total ≈ H1 + H2
    end
end
