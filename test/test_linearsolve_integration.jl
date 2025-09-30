"""
Test suite for LinearSolve.jl integration in ConvergentBornSolver.

This test suite validates that the LinearSolver object interface works correctly
and that different LinearSolver algorithms can be used.
"""

using Test
using FrequencyMaxwell
using LinearAlgebra
using LinearSolve: KrylovJL_GMRES, KrylovJL_BICGSTAB

"""
Test solver creation for validation tests with small problem sizes for speed.
"""
function create_test_solver(T = Float64; linear_solver = KrylovJL_GMRES())
    return ConvergentBornSolver(
        wavelength = T(500e-9),           # 500 nm
        permittivity_bg = T(1.33^2),      # Water background
        resolution = (T(50e-9), T(50e-9), T(50e-9)),  # 50 nm isotropic
        grid_size = (32, 32, 16),         # Small grid for fast testing
        boundary_thickness = (T(0.0), T(0.0), T(200e-9)),
        field_attenuation = (T(0.0), T(0.0), T(200e-9)),
        field_attenuation_sharpness = T(1.0),
        periodic_boundary = (true, true, false),
        iterations_max = 5,               # Small number for fast testing
        tolerance = T(1e-6),              # Reasonable tolerance for testing
        linear_solver = linear_solver
    )
end

"""
Create test permittivity distribution with known scattering features.
"""
function create_test_permittivity(solver::ConvergentBornSolver{T}) where {T}
    grid_size = solver.grid_size
    eps_bg = solver.permittivity_bg
    permittivity = fill(Complex{T}(eps_bg), grid_size)

    # Add a spherical scatterer in the center
    center = grid_size .÷ 2 .+ 1
    radius_pixels = 3  # Small scatterer

    for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
        dist_sq = (i - center[1])^2 + (j - center[2])^2 + (k - center[3])^2
        if dist_sq ≤ radius_pixels^2
            permittivity[i, j, k] = Complex{T}(1.46^2)  # SiO2 scatterer
        end
    end

    return permittivity
end

@testset "LinearSolve.jl Integration Tests" begin
    @testset "Configuration with LinearSolver Objects" begin
        @testset "GMRES Configuration" begin
            solver = create_test_solver(; linear_solver = KrylovJL_GMRES())
            @test isa(solver.linear_solver, typeof(KrylovJL_GMRES()))
            @test solver.wavelength ≈ 500e-9
        end

        @testset "BiCGSTAB Configuration" begin
            solver = create_test_solver(; linear_solver = KrylovJL_BICGSTAB())
            @test isa(solver.linear_solver, typeof(KrylovJL_BICGSTAB()))
        end

        @testset "Type Promotion" begin
            solver32 = create_test_solver(Float32; linear_solver = KrylovJL_GMRES())
            @test solver32.wavelength isa Float32
            @test solver32.tolerance isa Float32
            @test isa(solver32.linear_solver, typeof(KrylovJL_GMRES()))
        end
    end

    @testset "Solver Construction" begin
        solver_gmres = create_test_solver(; linear_solver = KrylovJL_GMRES())
        solver_bicgstab = create_test_solver(; linear_solver = KrylovJL_BICGSTAB())

        @test isa(solver_gmres.linear_solver, typeof(KrylovJL_GMRES()))
        @test isa(solver_bicgstab.linear_solver, typeof(KrylovJL_BICGSTAB()))
    end

    @testset "Basic Solver Functionality" begin
        @testset "GMRES Solver" begin
            solver = create_test_solver(; linear_solver = KrylovJL_GMRES())
            permittivity = create_test_permittivity(solver)

            # Create simple plane wave source
            source = PlaneWaveSource(
                wavelength = solver.wavelength,
                polarization = [1.0, 0.0, 0.0],
                k_vector = [0.0, 0.0, 1.0],
                amplitude = 1.0
            )

            # Test that solve doesn't error (basic functionality test)
            @test_nowarn E_field, H_field = solve(solver, source, permittivity)
        end

        @testset "BiCGSTAB Solver" begin
            solver = create_test_solver(; linear_solver = KrylovJL_BICGSTAB())
            permittivity = create_test_permittivity(solver)

            # Create simple plane wave source
            source = PlaneWaveSource(
                wavelength = solver.wavelength,
                polarization = [1.0, 0.0, 0.0],
                k_vector = [0.0, 0.0, 1.0],
                amplitude = 1.0
            )

            # Test that solve doesn't error (basic functionality test)
            @test_nowarn E_field, H_field = solve(solver, source, permittivity)
        end
    end

    @testset "Different LinearSolver Algorithms" begin
        algorithms = [KrylovJL_GMRES(), KrylovJL_BICGSTAB()]

        for algorithm in algorithms
            @testset "Algorithm: $(typeof(algorithm))" begin
                solver = create_test_solver(; linear_solver = algorithm)

                @test isa(solver.linear_solver, typeof(algorithm))

                # Test that configuration is valid
                @test solver.wavelength > 0
                @test solver.permittivity_bg > 0
                @test all(solver.resolution .> 0)
            end
        end
    end

    @testset "Field Output Validation" begin
        solver = create_test_solver(; linear_solver = KrylovJL_GMRES())
        permittivity = create_test_permittivity(solver)

        source = PlaneWaveSource(
            wavelength = solver.wavelength,
            polarization = [1.0, 0.0, 0.0],
            k_vector = [0.0, 0.0, 1.0],
            amplitude = 1.0
        )

        E_field, H_field = solve(solver, source, permittivity)

        # Test field properties
        @test size(E_field) == (solver.grid_size..., 3)
        @test size(H_field) == (solver.grid_size..., 3)
        @test eltype(E_field) <: Complex
        @test eltype(H_field) <: Complex

        # Test that fields are not trivial (should have some variation due to scatterer)
        @test !all(E_field .≈ E_field[1, 1, 1, 1])
        @test !all(H_field .≈ H_field[1, 1, 1, 1])

        # Test that fields are finite
        @test all(isfinite.(E_field))
        @test all(isfinite.(H_field))
    end
end
