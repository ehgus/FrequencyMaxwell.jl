"""
Basic solver functionality tests - migrated from performance_test.jl example
Tests fundamental solver initialization, configuration, and basic operations.
"""

@testset "Basic Solver Functionality" begin
    @testset "ConvergentBorn Solver Construction" begin
        # Test streamlined solver construction with keyword arguments
        solver = ConvergentBornSolver(
            wavelength = 532e-9,
            permittivity_bg = 1.333^2,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (64, 64, 32)
        )

        # Test solver properties
        @test typeof(solver) == ConvergentBornSolver{Float64}
        @test isnothing(solver.permittivity)  # Not set yet
        @test solver.Bornmax == 0  # Initial value

        # Test configuration fields (now directly accessible)
        @test solver.wavelength ≈ 532e-9
        @test solver.permittivity_bg ≈ 1.333^2
        @test solver.resolution == (50e-9, 50e-9, 50e-9)
        @test solver.grid_size == (64, 64, 32)

        # Test default values
        @test solver.use_abbe_sine == true
        @test solver.periodic_boundary == (true, true, false)
        @test solver.tolerance ≈ 1e-6
        @test solver.iterations_max == -1
    end

    @testset "Configuration Field Access" begin
        # Test direct field access and utility functions
        solver = ConvergentBornSolver(
            wavelength = 500e-9,
            permittivity_bg = 1.5^2,
            resolution = (25e-9, 25e-9, 25e-9),
            grid_size = (32, 32, 16),
            tolerance = 1e-4,
            iterations_max = 50
        )

        # Test direct field access
        @test solver.wavelength ≈ 500e-9
        @test solver.permittivity_bg ≈ 1.5^2
        @test solver.resolution == (25e-9, 25e-9, 25e-9)
        @test solver.grid_size == (32, 32, 16)
        @test solver.tolerance ≈ 1e-4
        @test solver.iterations_max == 50

        # Test utility functions
        @test grid_spacing(solver) == solver.resolution
        @test domain_size(solver) == (32 * 25e-9, 32 * 25e-9, 16 * 25e-9)
        @test wavenumber_background(solver) ≈ 2π * sqrt(1.5^2) / 500e-9
    end

    @testset "Type Promotion" begin
        # Test that different numeric types are properly promoted
        solver = ConvergentBornSolver(
            wavelength = 532,  # Int -> Float64
            permittivity_bg = 1.333f0^2,  # Float32
            resolution = (50.0, 50.0, 50.0),  # Float64
            grid_size = (64, 64, 32)
        )

        @test typeof(solver.wavelength) == Float64
        @test typeof(solver.permittivity_bg) == Float64
        @test eltype(solver.resolution) == Float64
    end

    @testset "Solver State Management" begin
        solver = ConvergentBornSolver(
            wavelength = 500e-9,
            permittivity_bg = 1.33^2,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 16)
        )

        # Test initial state
        @test solver.iteration_count == 0
        @test isempty(solver.residual_history)
        @test solver.Bornmax == 0

        # Test reset functionality
        reset!(solver)
        @test solver.iteration_count == 0
        @test isempty(solver.residual_history)
    end

    @testset "Solver Parameter Validation" begin
        # Test parameter validation during construction
        @test_throws ArgumentError ConvergentBornSolver(
            wavelength = -500e-9,  # Invalid: negative wavelength
            permittivity_bg = 1.0,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 32)
        )

        @test_throws ArgumentError ConvergentBornSolver(
            wavelength = 500e-9,
            permittivity_bg = -1.0,  # Invalid: negative permittivity
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 32)
        )

        @test_throws ArgumentError ConvergentBornSolver(
            wavelength = 500e-9,
            permittivity_bg = 1.0,
            resolution = (-50e-9, 50e-9, 50e-9),  # Invalid: negative resolution
            grid_size = (32, 32, 32)
        )

        # Test valid construction should work
        solver = ConvergentBornSolver(
            wavelength = 633e-9,
            permittivity_bg = 1.0,  # Vacuum
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 32)
        )

        # Test that solver was created successfully
        @test solver.permittivity_bg ≈ 1.0
        @test solver.wavelength ≈ 633e-9
    end
end
