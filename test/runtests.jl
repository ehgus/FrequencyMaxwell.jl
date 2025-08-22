using Test
using FrequencyMaxwell

@testset "FrequencyMaxwell.jl Tests" begin
    @testset "Core Configuration" begin
        @testset "ConvergentBornConfig Construction" begin
            # Test basic construction
            config = ConvergentBornConfig(
                wavelength = 500e-9,
                NA = 1.4,
                permittivity_bg = 1.33^2,
                resolution = (50e-9, 50e-9, 50e-9),
                grid_size = (128, 128, 32)
            )
            
            @test config.wavelength ≈ 500e-9
            @test config.NA ≈ 1.4
            @test config.permittivity_bg ≈ 1.33^2
            @test config.resolution == (50e-9, 50e-9, 50e-9)
            @test config.grid_size == (128, 128, 32)
            
            # Test validation
            @test_throws ArgumentError ConvergentBornConfig(
                wavelength = -1.0,  # Invalid
                NA = 1.4,
                resolution = (50e-9, 50e-9, 50e-9),
                grid_size = (128, 128, 32)
            )
            
            @test_throws ArgumentError ConvergentBornConfig(
                wavelength = 500e-9,
                NA = 3.0,  # Invalid (> 2)
                resolution = (50e-9, 50e-9, 50e-9),
                grid_size = (128, 128, 32)
            )
        end
        
        @testset "Configuration Utilities" begin
            config = ConvergentBornConfig(
                wavelength = 500e-9,
                NA = 1.4,
                resolution = (50e-9, 50e-9, 100e-9),
                grid_size = (64, 64, 32)
            )
            
            # Test domain size calculation
            domain = domain_size(config)
            @test domain[1] ≈ 64 * 50e-9
            @test domain[2] ≈ 64 * 50e-9  
            @test domain[3] ≈ 32 * 100e-9
            
            # Test wavenumber calculation
            k0 = wavenumber_background(config)
            @test k0 ≈ 2π * sqrt(config.permittivity_bg) / config.wavelength
        end
    end
    
    @testset "Solver Construction" begin
        config = ConvergentBornConfig(
            wavelength = 500e-9,
            NA = 1.4,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 16)
        )
        
        solver = ConvergentBornSolver(config)
        
        @test solver.config === config
        @test solver.iteration_count == 0
        @test isempty(solver.residual_history)
        @test solver.Green_function === nothing
        @test solver.potential === nothing
    end
    
    @testset "Source Types" begin
        @testset "PlaneWaveSource Construction" begin
            # Test basic construction
            source = PlaneWaveSource(
                wavelength = 500e-9,
                polarization = [1.0, 0.0, 0.0],
                k_vector = [0.0, 0.0, 1.0],
                amplitude = 1.0
            )
            
            @test source.wavelength ≈ 500e-9
            @test source.amplitude ≈ 1.0
            @test source.phase ≈ 0.0
            
            # Test validation - transverse condition
            @test_throws ArgumentError PlaneWaveSource(
                wavelength = 500e-9,
                polarization = [1.0, 0.0, 1.0],  # Not transverse
                k_vector = [0.0, 0.0, 1.0],
                amplitude = 1.0
            )
            
            # Test validation - normalized k-vector
            @test_throws ArgumentError PlaneWaveSource(
                wavelength = 500e-9,
                polarization = [1.0, 0.0, 0.0],
                k_vector = [0.0, 0.0, 2.0],  # Not normalized
                amplitude = 1.0
            )
        end
        
        @testset "Source Interface Methods" begin
            source = PlaneWaveSource(
                wavelength = 633e-9,
                polarization = [1.0, 0.0, 0.0],
                k_vector = [0.0, 0.0, 1.0],
                amplitude = 2.0
            )
            
            # Test interface methods
            @test source_wavelength(source) ≈ 633e-9
            @test source_power(source) > 0
            @test validate_source(source) == true
        end
    end
    
    @testset "Phantom Generation" begin
        @testset "Bead Phantom" begin
            phantom = phantom_bead(
                (64, 64, 32),           # Grid size
                [1.5^2],                # Permittivity
                4.0                     # Radius
            )
            
            @test size(phantom) == (64, 64, 32)
            @test eltype(phantom) <: Complex
            
            # Check that bead has higher permittivity than background
            center_idx = (32, 32, 16)
            @test real(phantom[center_idx...]) > 1.0
            
            # Check background
            corner_idx = (1, 1, 1)
            @test real(phantom[corner_idx...]) ≈ 1.0
        end
        
        @testset "Plate Phantom" begin
            phantom = phantom_plate(
                (64, 64, 32),           # Grid size
                [1.8^2],                # Permittivity
                8.0                     # Thickness
            )
            
            @test size(phantom) == (64, 64, 32)
            @test eltype(phantom) <: Complex
            
            # Check plate region (Z-normal by default)
            center_z = 16
            @test real(phantom[32, 32, center_z]) > 1.0
            
            # Check outside plate
            @test real(phantom[32, 32, 1]) ≈ 1.0
            @test real(phantom[32, 32, 32]) ≈ 1.0
        end
    end
    
    @testset "Electromagnetic Fields" begin
        # Create test field arrays
        grid_size = (32, 32, 16)
        E_field = zeros(ComplexF64, grid_size..., 3)
        H_field = zeros(ComplexF64, grid_size..., 3)
        resolution = (50e-9, 50e-9, 100e-9)
        wavelength = 500e-9
        
        # Add some test data
        E_field[:, :, :, 1] .= 1.0  # Ex = 1
        H_field[:, :, :, 2] .= 1.0/377.0  # Hy = Ex/Z0
        
        fields = ElectromagneticField(E_field, H_field, grid_size, resolution, wavelength)
        
        @test fields.grid_size == grid_size
        @test fields.resolution == resolution
        @test fields.wavelength == wavelength
        
        # Test derived quantities
        domain = domain_size(fields)
        @test domain[1] ≈ 32 * 50e-9
        @test domain[2] ≈ 32 * 50e-9
        @test domain[3] ≈ 16 * 100e-9
        
        # Test field energy calculation
        energy = field_energy(fields)
        @test energy > 0
        
        # Test intensity calculation
        intensity = field_intensity(fields)
        @test size(intensity) == grid_size
        @test all(intensity .≈ 1.0)  # |Ex|² = 1
    end
    
    @testset "Integration Tests" begin
        @testset "Basic Solver Pipeline" begin
            # Create configuration
            config = ConvergentBornConfig(
                wavelength = 500e-9,
                NA = 1.4,
                resolution = (100e-9, 100e-9, 100e-9),
                grid_size = (16, 16, 8)  # Small for fast testing
            )
            
            # Create solver
            solver = ConvergentBornSolver(config)
            
            # Create source
            source = PlaneWaveSource(
                wavelength = 500e-9,
                polarization = [1.0, 0.0, 0.0],
                k_vector = [0.0, 0.0, 1.0]
            )
            
            # Create simple phantom
            permittivity = ones(ComplexF64, config.grid_size...)
            permittivity[8, 8, 4] = 1.5^2  # Single high-index voxel
            
            # Test that solve runs without error (placeholder implementation)
            @test_nowarn begin
                E_field, H_field = solve(solver, source, permittivity)
                @test size(E_field) == (config.grid_size..., 3)
                @test size(H_field) == (config.grid_size..., 3)
            end
        end
    end
end
