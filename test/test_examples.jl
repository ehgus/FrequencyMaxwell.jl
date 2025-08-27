"""
Tests for FrequencyMaxwell examples to ensure they run without error.
This migrates the working examples into the standardized test suite format.
"""

using Test
using FrequencyMaxwell

@testset "Example Integration Tests" begin
    
    @testset "Basic Configuration Examples" begin
        # Test basic scattering configuration (from basic_scattering.jl example)
        @test_nowarn begin
            config = ConvergentBornConfig(
                wavelength = 500e-9,      # 500 nm (green light)
                permittivity_bg = 1.33^2, # Water background (n=1.33)
                resolution = (50e-9, 50e-9, 50e-9),  # 50 nm isotropic resolution
                grid_size = (64, 64, 32)  # Smaller grid for faster testing
            )
            
            @test config.wavelength == 500e-9
            @test config.permittivity_bg == 1.33^2
        end
        
        # Test solver creation (basic example workflow)
        @test_nowarn begin
            config = ConvergentBornConfig(
                wavelength = 632.8e-9,   # He-Ne laser (from examples)
                permittivity_bg = 1.0,   # Air background
                resolution = (100e-9, 100e-9, 100e-9),
                grid_size = (32, 32, 16)
            )
            
            solver = ConvergentBornSolver(config)
            @test solver isa ConvergentBornSolver
        end
    end
    
    @testset "Phantom Generation Examples" begin
        # Test phantom generation (from phantom_gallery.jl example)
        @test_nowarn begin
            grid_size = (32, 32, 16)  # Small for testing
            
            # Test bead phantom (key example from gallery)
            bead_phantom = phantom_bead(grid_size, [2.25], 5)
            @test size(bead_phantom) == grid_size
            @test eltype(bead_phantom) == ComplexF64
        end
        
        @test_nowarn begin
            grid_size = (32, 32, 16)
            
            # Test cylinder phantom (key example from gallery) 
            cyl_phantom = phantom_cylinder(grid_size, [3.0], 4)
            @test size(cyl_phantom) == grid_size
            @test eltype(cyl_phantom) == ComplexF64
        end
        
        @test_nowarn begin
            grid_size = (32, 32, 16)
            
            # Test plate phantom (key example from gallery)
            plate_phantom = phantom_plate(grid_size, [2.25], 8)
            @test size(plate_phantom) == grid_size
            @test eltype(plate_phantom) == ComplexF64
        end
    end
    
    @testset "Field Source Examples" begin
        # Test plane wave source (used in most examples) - with correct API
        @test_nowarn begin
            config = ConvergentBornConfig(
                wavelength = 500e-9,
                permittivity_bg = 1.0,
                resolution = (100e-9, 100e-9, 100e-9),
                grid_size = (16, 16, 8)
            )
            
            # Test plane wave source with correct keyword arguments
            source = PlaneWaveSource(
                wavelength = config.wavelength,
                polarization = [1.0, 0.0, 0.0],  # x-polarized
                k_vector = [0.0, 0.0, 1.0], # z-propagating
                amplitude = 1.0,
                phase = 0.0
            )
            
            @test source isa PlaneWaveSource
            @test source_wavelength(source) ≈ config.wavelength
        end
        
        # Test incident field generation (part of solve workflow in examples)
        @test_nowarn begin
            config = ConvergentBornConfig(
                wavelength = 532e-9,
                permittivity_bg = 1.0,
                resolution = (200e-9, 200e-9, 200e-9),
                grid_size = (16, 16, 8)
            )
            
            source = PlaneWaveSource(
                wavelength = config.wavelength,
                polarization = [1.0, 0.0, 0.0],  # x-polarized
                k_vector = [0.0, 0.0, 1.0], # z-propagating
                amplitude = 1.0,
                phase = 0.0
            )
            
            incident_fields = generate_incident_fields(source, config)
            @test incident_fields isa ElectromagneticField
            @test size(incident_fields.E) == (16, 16, 8, 3)
            @test size(incident_fields.H) == (16, 16, 8, 3)
        end
    end
    
    @testset "Configuration Utilities Examples" begin
        # Test utility functions used in examples
        @test_nowarn begin
            config = ConvergentBornConfig(
                wavelength = 532e-9,
                resolution = (50e-9, 50e-9, 100e-9),
                grid_size = (32, 32, 16)
            )
            
            # Domain size calculation (used in many examples)
            domain = domain_size(config)
            @test length(domain) == 3
            @test domain[1] ≈ 32 * 50e-9
            @test domain[2] ≈ 32 * 50e-9
            @test domain[3] ≈ 16 * 100e-9
            
            # Grid spacing (utility used in examples)
            spacing = grid_spacing(config)
            @test spacing == config.resolution
            
            # Background wavenumber (physics calculation in examples)
            k_bg = wavenumber_background(config)
            @test k_bg ≈ 2π * sqrt(config.permittivity_bg) / config.wavelength
        end
    end
    
    @testset "Field Analysis Examples" begin
        # Test field analysis utilities (used in examples for post-processing)
        @test_nowarn begin
            # Create test electromagnetic field (typical example result)
            E_data = ones(ComplexF64, 16, 16, 8, 3) * 0.1
            H_data = ones(ComplexF64, 16, 16, 8, 3) * 0.01
            
            fields = ElectromagneticField(
                E_data, H_data, (16, 16, 8), (100e-9, 100e-9, 200e-9), 633e-9
            )
            
            # Field energy calculation (common in examples)
            energy = field_energy(fields)
            @test energy > 0
            @test isfinite(energy)
            
            # Intensity calculation (visualization in examples)
            intensity = field_intensity(fields)
            @test size(intensity) == (16, 16, 8)
            @test all(intensity .≥ 0)
            
            # Poynting vector (power flow analysis in examples) - returns complex
            S = poynting_vector(fields)
            @test size(S) == (16, 16, 8, 3)
            @test eltype(S) <: Complex  # Fixed: should be Complex, not Real
        end
    end
    
    @testset "Example Error Validation" begin
        # Test error handling that examples should implement
        
        # Invalid wavelength (example validation)
        @test_throws ArgumentError ConvergentBornConfig(
            wavelength = -500e-9,  # Invalid
            permittivity_bg = 1.33^2,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 16)
        )
    end
    
    @testset "Example Integration Workflow" begin
        # Test a complete example workflow (simplified from actual examples)
        @test_nowarn begin
            # Step 1: Configuration (from examples)
            config = ConvergentBornConfig(
                wavelength = 532e-9,  # Green laser (common in examples)
                permittivity_bg = 1.333^2,  # Water (common background)
                resolution = (100e-9, 100e-9, 100e-9),
                grid_size = (24, 24, 12)  # Small for testing
            )
            
            # Step 2: Phantom creation (from phantom_gallery.jl)
            phantom = phantom_bead(config.grid_size, [1.4607^2], 3)  # SiO2 bead
            
            # Step 3: Source creation (from examples) - with correct API
            source = PlaneWaveSource(
                wavelength = config.wavelength,
                polarization = [1.0, 0.0, 0.0],  # x-polarized
                k_vector = [0.0, 0.0, 1.0], # z-propagating
                amplitude = 1.0,
                phase = 0.0
            )
            
            # Step 4: Solver creation (common workflow)
            solver = ConvergentBornSolver(config)
            
            # Verify all components were created successfully
            @test size(phantom) == config.grid_size
            @test source isa PlaneWaveSource
            @test solver isa ConvergentBornSolver
        end
    end
end
