"""
Field source tests - migrated from forward solver and grating examples
Tests electromagnetic source generation and field initialization.
"""

@testset "Electromagnetic Sources" begin
    @testset "Plane Wave Source" begin
        # Test basic plane wave source - matches original examples
        config = ConvergentBornConfig(
            wavelength = 532e-9,
            NA = 1.2,
            permittivity_bg = 1.333^2,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (64, 64, 32)
        )
        
        # Test plane wave propagating in z-direction
        source = PlaneWaveSource(
            config,
            polarization = (1.0, 0.0, 0.0),  # x-polarized
            direction = 3,  # z-direction
            k_transverse = (0.0, 0.0)  # Normal incidence
        )
        
        @test source.config == config
        @test source.polarization == (1.0, 0.0, 0.0)
        @test source.direction == 3
        @test source.k_transverse == (0.0, 0.0)
        
        # Test field generation
        E_field, H_field = generate_fields(source)
        
        @test size(E_field) == (config.grid_size..., 3)
        @test size(H_field) == (config.grid_size..., 3)
        
        # Test field properties
        @test !any(isnan.(E_field))
        @test !any(isnan.(H_field))
        @test !any(isinf.(E_field))
        @test !any(isinf.(H_field))
    end
    
    @testset "Source Polarization" begin
        config = ConvergentBornConfig(
            wavelength = 500e-9,
            NA = 1.0,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 16)
        )
        
        # Test x-polarization
        source_x = PlaneWaveSource(
            config,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )
        E_x, H_x = generate_fields(source_x)
        
        # Test y-polarization
        source_y = PlaneWaveSource(
            config,
            polarization = (0.0, 1.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )
        E_y, H_y = generate_fields(source_y)
        
        # Test circular polarization
        source_circ = PlaneWaveSource(
            config,
            polarization = (1.0, 1.0im, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )
        E_circ, H_circ = generate_fields(source_circ)
        
        # Verify orthogonality for linear polarizations
        center = div.(config.grid_size, 2) .+ 1
        Ex_center = E_x[center..., 1]
        Ey_center = E_y[center..., 2]
        @test abs(real(Ex_center * conj(Ey_center))) < 1e-10  # Should be orthogonal
    end
    
    @testset "Oblique Incidence" begin
        config = ConvergentBornConfig(
            wavelength = 355e-9,  # UV wavelength from grating example
            NA = 1.0,
            permittivity_bg = 1.0,
            resolution = (10e-9, 10e-9, 10e-9),
            grid_size = (20, 50, 50)
        )
        
        # Test different diffraction orders - based on grating example
        illumination_orders = [-3, -1, 0, 1, 3]
        
        for order in illumination_orders
            # Calculate transverse k-vector for this order
            k_y = 2π * order / (config.grid_size[2] * config.resolution[2])
            k_transverse = (0.0, k_y)
            
            source = PlaneWaveSource(
                config,
                polarization = (1.0, 0.0, 0.0),
                direction = 3,
                k_transverse = k_transverse
            )
            
            E_field, H_field = generate_fields(source)
            
            # Test that fields are generated properly
            @test size(E_field) == (config.grid_size..., 3)
            @test size(H_field) == (config.grid_size..., 3)
            
            # Test that oblique incidence creates expected phase patterns
            if order != 0
                # Should have phase variation across y-direction
                phase_y1 = angle(E_field[10, 1, 25, 1])
                phase_y2 = angle(E_field[10, end, 25, 1])
                @test abs(phase_y1 - phase_y2) > 0.1  # Significant phase difference
            end
        end
    end
    
    @testset "Source Field Properties" begin
        config = ConvergentBornConfig(
            wavelength = 633e-9,
            NA = 1.2,
            permittivity_bg = 1.33^2,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 32)
        )
        
        source = PlaneWaveSource(
            config,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )
        
        E_field, H_field = generate_fields(source)
        
        # Test Maxwell's equations in source region
        # For plane wave: ∇ × E = -iωμH, ∇ × H = iωεE
        
        # Test field magnitudes are reasonable
        E_magnitude = sqrt.(sum(abs2.(E_field), dims=4))
        H_magnitude = sqrt.(sum(abs2.(H_field), dims=4))
        
        @test maximum(E_magnitude) > 0
        @test maximum(H_magnitude) > 0
        
        # Test impedance relationship for plane wave
        # |H| ≈ |E| * sqrt(ε/μ) = |E| * sqrt(ε₀/μ₀) * sqrt(εᵣ)
        n_bg = sqrt(config.permittivity_bg)
        Z0 = 377.0  # Free space impedance
        Z_medium = Z0 / n_bg
        
        center = div.(config.grid_size, 2) .+ 1
        E_center_mag = abs(E_field[center..., 1])
        H_center_mag = abs(H_field[center..., 2])  # H_y for E_x
        
        if H_center_mag > 0
            impedance_ratio = E_center_mag / H_center_mag
            @test impedance_ratio ≈ Z_medium rtol=0.1
        end
    end
    
    @testset "Multiple Sources" begin
        config = ConvergentBornConfig(
            wavelength = 500e-9,
            NA = 1.0,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32)
        )
        
        # Create multiple sources for coherent superposition
        sources = [
            PlaneWaveSource(config, polarization=(1.0, 0.0, 0.0), direction=3, k_transverse=(0.0, 0.0)),
            PlaneWaveSource(config, polarization=(0.0, 1.0, 0.0), direction=3, k_transverse=(0.0, 0.0))
        ]
        
        E_total, H_total = superpose_sources(sources)
        
        @test size(E_total) == (config.grid_size..., 3)
        @test size(H_total) == (config.grid_size..., 3)
        
        # Test that superposition is linear
        E1, H1 = generate_fields(sources[1])
        E2, H2 = generate_fields(sources[2])
        
        @test E_total ≈ E1 + E2
        @test H_total ≈ H1 + H2
    end
end
