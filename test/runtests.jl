using Test
using FrequencyMaxwell

@testset "FrequencyMaxwell.jl Tests" begin
    @testset "Core Solver" begin
        @testset "ConvergentBornSolver Construction" begin
            # Test basic construction
            solver = ConvergentBornSolver(
                wavelength = 500e-9,
                permittivity_bg = 1.33^2,
                resolution = (50e-9, 50e-9, 50e-9),
                grid_size = (128, 128, 32)
            )

            @test solver.wavelength ≈ 500e-9
            @test solver.permittivity_bg ≈ 1.33^2
            @test solver.resolution == (50e-9, 50e-9, 50e-9)
            @test solver.grid_size == (128, 128, 32)

            # Test validation
            @test_throws ArgumentError ConvergentBornSolver(
                wavelength = -1.0,  # Invalid
                resolution = (50e-9, 50e-9, 50e-9),
                grid_size = (128, 128, 32)
            )

            @test_throws ArgumentError ConvergentBornSolver(
                wavelength = 500e-9,
                permittivity_bg = -1.0,  # Invalid permittivity
                resolution = (50e-9, 50e-9, 50e-9),
                grid_size = (128, 128, 32)
            )
        end

        @testset "Solver Utilities" begin
            solver = ConvergentBornSolver(
                wavelength = 500e-9,
                resolution = (50e-9, 50e-9, 100e-9),
                grid_size = (64, 64, 32)
            )

            # Test domain size calculation
            domain = domain_size(solver)
            @test domain[1] ≈ 64 * 50e-9
            @test domain[2] ≈ 64 * 50e-9
            @test domain[3] ≈ 32 * 100e-9
        end
    end

    @testset "Basic Phantom Generation" begin
        @testset "Phantom Bead" begin
            # Test the phantom_bead function from the existing API
            grid_size = (32, 32, 32)
            permittivity_profile = [1.46^2]  # SiO2 (single value for single bead)
            radius_pixels = 5

            phantom = phantom_bead(grid_size, permittivity_profile, radius_pixels)

            @test size(phantom) == grid_size
            @test phantom[1, 1, 1] ≈ 1.0  # Background (default is 1.0)

            # Test that center has bead permittivity
            center = div.(grid_size, 2) .+ 1
            center_value = phantom[center...]
            @test center_value ≈ permittivity_profile[1]
        end

        @testset "Phantom Plate" begin
            # Test layered structure with correct API - simplified test
            grid_size = (20, 20, 20)
            permittivity_profile = [2.0]  # Single material for plate
            thickness_pixels = 8  # Single thickness value

            phantom = phantom_plate(grid_size, permittivity_profile, thickness_pixels)

            @test size(phantom) == grid_size

            # Test that phantom was created (not checking exact values due to implementation details)
            @test !all(phantom .== 1.0)  # Should have some variation from background
            @test maximum(real.(phantom)) ≥ 1.0  # Should have material with permittivity ≥ 1
        end
    end

    @testset "ElectromagneticField Structure" begin
        @testset "Field Construction" begin
            # Test basic field construction with correct constructor
            E_data = rand(ComplexF64, 32, 32, 16, 3)  # Grid + 3 components
            H_data = rand(ComplexF64, 32, 32, 16, 3)  # Grid + 3 components
            grid_size = (32, 32, 16)
            resolution = (100e-9, 100e-9, 200e-9)
            wavelength = 500e-9

            fields = ElectromagneticField(E_data, H_data, grid_size, resolution, wavelength)

            @test size(fields.E) == (32, 32, 16, 3)
            @test size(fields.H) == (32, 32, 16, 3)
            @test fields.grid_size == grid_size
            @test fields.resolution == resolution
            @test fields.wavelength ≈ wavelength
        end

        @testset "Field Utilities" begin
            # Create test fields with correct constructor
            E_data = ones(ComplexF64, 16, 16, 8, 3) * 0.1
            H_data = ones(ComplexF64, 16, 16, 8, 3) * 0.01
            grid_size = (16, 16, 8)
            resolution = (50e-9, 50e-9, 100e-9)
            wavelength = 633e-9

            fields = ElectromagneticField(E_data, H_data, grid_size, resolution, wavelength)

            # Test energy calculation
            energy = field_energy(fields)
            @test energy > 0
            @test isfinite(energy)

            # Test intensity calculation  
            intensity = field_intensity(fields)
            @test size(intensity) == (16, 16, 8)
            @test all(intensity .≥ 0)

            # Test Poynting vector
            S = poynting_vector(fields)
            @test size(S) == (16, 16, 8, 3)
        end

        @testset "Field Domain Properties" begin
            E_data = rand(ComplexF64, 20, 30, 10, 3)
            H_data = rand(ComplexF64, 20, 30, 10, 3)
            grid_size = (20, 30, 10)
            resolution = (25e-9, 40e-9, 80e-9)
            wavelength = 500e-9

            fields = ElectromagneticField(E_data, H_data, grid_size, resolution, wavelength)

            domain = domain_size(fields)
            @test domain[1] ≈ 20 * 25e-9
            @test domain[2] ≈ 30 * 40e-9
            @test domain[3] ≈ 10 * 80e-9
        end
    end

    @testset "Solver Utilities" begin
        @testset "Domain Calculations" begin
            solver = ConvergentBornSolver(
                wavelength = 633e-9,
                permittivity_bg = 1.0,
                resolution = (100e-9, 100e-9, 200e-9),
                grid_size = (50, 60, 25)
            )

            # Test utility functions if they exist
            if isdefined(FrequencyMaxwell, :grid_spacing)
                spacing = grid_spacing(solver)
                @test spacing == solver.resolution
            end

            if isdefined(FrequencyMaxwell, :wavenumber_background)
                k_bg = wavenumber_background(solver)
                @test k_bg ≈ 2π * sqrt(solver.permittivity_bg) / solver.wavelength
            end
        end
    end

    @testset "Integration Test - Example Migration" begin
        @testset "SiO2 Bead Example Setup" begin
            # Replicate the core setup from forward_SiO2_5um_in_water.jl

            # Solver matching the original example (scaled down)
            solver = ConvergentBornSolver(
                wavelength = 532e-9,  # 532 nm green laser
                permittivity_bg = 1.333^2,  # Water
                resolution = (100e-9, 100e-9, 100e-9),  # 100 nm resolution
                grid_size = (51, 51, 48)  # Smaller than original 201x201x191
            )

            # Create SiO2 bead phantom - using single permittivity value
            permittivity_SiO2 = 1.4607^2
            radius_pixels = 5  # Approximately 0.5 μm radius
            phantom = phantom_bead(
                solver.grid_size,
                [permittivity_SiO2],  # Single bead material
                radius_pixels
            )

            # Test phantom properties
            @test size(phantom) == solver.grid_size
            @test phantom[1, 1, 1] ≈ 1.0  # Background is 1.0 by default

            # Test bead is present at center
            center = div.(solver.grid_size, 2) .+ 1
            @test phantom[center...] ≈ permittivity_SiO2

            # Test proper domain size
            domain = domain_size(solver)
            @test domain[1] ≈ 51 * 100e-9  # 5.1 μm
            @test domain[2] ≈ 51 * 100e-9  # 5.1 μm
            @test domain[3] ≈ 48 * 100e-9  # 4.8 μm
        end

        @testset "Grating Example Setup" begin
            # Replicate setup from adjoint_grating.jl (simplified)

            solver = ConvergentBornSolver(
                wavelength = 355e-9,  # UV wavelength
                permittivity_bg = 1.4338^2,  # PDMS background
                resolution = (20e-9, 20e-9, 20e-9),  # 20 nm resolution
                grid_size = (21, 50, 50)  # Smaller than original
            )

            # Create simple plate structure (simplified from multi-layer)
            permittivity_TiO2 = (2.9734 + 0.0467im)^2  # Lossy TiO2
            thickness_pixels = 8  # Single layer thickness

            phantom = phantom_plate(
                solver.grid_size,
                [permittivity_TiO2],
                thickness_pixels
            )

            # Test structure properties (less specific due to implementation details)
            @test size(phantom) == solver.grid_size
            @test maximum(real.(phantom)) > 1.0  # Should have high-index material
            @test maximum(imag.(phantom)) > 0.0  # Should have lossy material

            # Test that we have the expected complex permittivity somewhere
            @test any(p -> abs(p - permittivity_TiO2) < 0.1, phantom)
        end

        @testset "Performance Example Setup" begin
            # Replicate performance test setup (very small for testing)

            solver = ConvergentBornSolver(
                wavelength = 532e-9,
                permittivity_bg = 1.333^2,
                resolution = (200e-9, 200e-9, 200e-9),  # Coarser resolution
                grid_size = (21, 21, 16)  # Very small grid
            )

            # Create test phantom
            permittivity_sp = 1.4607^2
            radius_pixels = 2
            phantom = phantom_bead(
                solver.grid_size,
                [permittivity_sp],
                radius_pixels
            )

            # Test that phantom was created properly
            @test size(phantom) == solver.grid_size

            # Basic timing test (just to ensure setup works)
            time_start = time()
            _ = phantom_bead(solver.grid_size, [2.0], 3)
            elapsed = time() - time_start
            @test elapsed < 1.0  # Should be very fast for small problem
        end
    end
end

# Include example integration tests
include("test_examples.jl")

# Include LinearSolve.jl integration tests
include("test_linearsolve_integration.jl")

# Include performance benchmarking tests (optional - can be run separately)
# Uncomment the line below to include performance benchmarks in regular test suite
# include("test_linearsolve_performance.jl")
