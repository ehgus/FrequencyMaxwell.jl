"""
Basic solver functionality tests - migrated from performance_test.jl example
Tests fundamental solver initialization, configuration, and basic operations.
"""

@testset "Basic Solver Functionality" begin
    @testset "ConvergentBorn Solver Construction" begin
        # Test streamlined solver construction with keyword arguments
        solver = ConvergentBornSolver(
            permittivity_bg = 1.333^2,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (64, 64, 32),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        # Test solver properties
        @test typeof(solver) == ConvergentBornSolver{Float64}
        @test isnothing(solver.permittivity)  # Not set yet
        @test solver.Bornmax == 0  # Initial value

        # Test configuration fields (now directly accessible)
        @test solver.permittivity_bg ≈ 1.333^2
        @test solver.resolution == (50e-9, 50e-9, 50e-9)
        @test solver.grid_size == (64, 64, 32)

        # Test boundary conditions
        @test length(solver.boundary_conditions) == 3
        @test all(bc -> bc isa PeriodicBoundaryCondition, solver.boundary_conditions)

        # Test default values
        @test solver.tolerance ≈ 1e-6
        @test solver.iterations_max == -1
    end

    @testset "Configuration Field Access" begin
        # Test direct field access and utility functions
        solver = ConvergentBornSolver(
            permittivity_bg = 1.5^2,
            resolution = (25e-9, 25e-9, 25e-9),
            grid_size = (32, 32, 16),
            boundary_conditions = AbsorbingBoundaryCondition(thickness = 1e-6, attenuation_thickness = 0.8e-6),
            tolerance = 1e-4,
            iterations_max = 50
        )

        # Test direct field access
        @test solver.permittivity_bg ≈ 1.5^2
        @test solver.resolution == (25e-9, 25e-9, 25e-9)
        @test solver.grid_size == (32, 32, 16)
        @test solver.tolerance ≈ 1e-4
        @test solver.iterations_max == 50

        # Test utility functions
        @test grid_spacing(solver) == solver.resolution
        @test domain_size(solver) == (32 * 25e-9, 32 * 25e-9, 16 * 25e-9)
    end

    @testset "Type Promotion" begin
        # Test that different numeric types are properly promoted
        solver = ConvergentBornSolver(
            permittivity_bg = 1.333f0^2,  # Float32
            resolution = (50.0, 50.0, 50.0),  # Float64
            grid_size = (64, 64, 32),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        @test typeof(solver.permittivity_bg) == Float64
        @test eltype(solver.resolution) == Float64
    end

    @testset "Solver State Management" begin
        solver = ConvergentBornSolver(
            permittivity_bg = 1.33^2,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 16),
            boundary_conditions = PeriodicBoundaryCondition()
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
            permittivity_bg = -1.0,  # Invalid: negative permittivity
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 32),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        @test_throws ArgumentError ConvergentBornSolver(
            permittivity_bg = 1.0,
            resolution = (-50e-9, 50e-9, 50e-9),  # Invalid: negative resolution
            grid_size = (32, 32, 32),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        # Test valid construction should work
        solver = ConvergentBornSolver(
            permittivity_bg = 1.0,  # Vacuum
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 32),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        # Test that solver was created successfully
        @test solver.permittivity_bg ≈ 1.0
    end

    @testset "Domain Size Calculations" begin
        solver = ConvergentBornSolver(
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 200e-9),
            grid_size = (50, 60, 25),
            boundary_conditions = PeriodicBoundaryCondition()
        )

        # Test domain size calculation
        domain = domain_size(solver)
        @test domain[1] ≈ 50 * 100e-9
        @test domain[2] ≈ 60 * 100e-9
        @test domain[3] ≈ 25 * 200e-9

        # Test utility functions if they exist
        if isdefined(FrequencyMaxwell, :grid_spacing)
            spacing = grid_spacing(solver)
            @test spacing == solver.resolution
        end
    end
end
