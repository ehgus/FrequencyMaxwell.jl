"""
Basic solver functionality tests - migrated from performance_test.jl example
Tests fundamental solver initialization, configuration, and basic operations.
"""

@testset "Basic Solver Functionality" begin
    @testset "ConvergentBorn Solver Construction" begin
        # Test basic solver construction with minimal parameters
        config = ConvergentBornConfig(
            wavelength = 532e-9,  # 532 nm (green laser)
            permittivity_bg = 1.333^2,  # Water background
            resolution = (50e-9, 50e-9, 50e-9),  # 50 nm resolution
            grid_size = (64, 64, 32)  # Smaller grid for testing
        )

        solver = ConvergentBornSolver(config)

        @test solver.config == config
        @test isnothing(solver.permittivity)  # Not set yet
        @test solver.Born_max == 100  # Default value

        # Test configuration access
        @test solver.config.wavelength ≈ 532e-9
        @test solver.config.permittivity_bg ≈ 1.333^2
    end

    @testset "Solver State Management" begin
        config = ConvergentBornConfig(
            wavelength = 500e-9,
            permittivity_bg = 1.33^2,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 16)
        )

        solver = ConvergentBornSolver(config)

        # Test Born iteration limits
        set_Born_max!(solver, 50)
        @test solver.Born_max == 50

        # Test convergence tolerance if implemented
        if hasfield(typeof(solver), :convergence_tolerance)
            set_convergence_tolerance!(solver, 1e-6)
            @test solver.convergence_tolerance ≈ 1e-6
        end
    end

    @testset "Solver Parameter Validation" begin
        config = ConvergentBornConfig(
            wavelength = 633e-9,
            permittivity_bg = 1.0,  # Vacuum
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 32)
        )

        solver = ConvergentBornSolver(config)

        # Test that solver handles edge cases properly
        @test solver.config.permittivity_bg ≈ 1.0

        # Test Born max validation
        @test_throws ArgumentError set_Born_max!(solver, -1)
        @test_throws ArgumentError set_Born_max!(solver, 0)
    end
end
