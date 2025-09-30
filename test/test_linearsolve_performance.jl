"""
Performance benchmarking suite for LinearSolve.jl integration.

This test suite measures and validates the performance of different LinearSolver
algorithms with the current FrequencyMaxwell interface.
"""

using Test
using FrequencyMaxwell
using LinearAlgebra

"""
Performance test solver creation with medium problem sizes for benchmarking.
"""
function create_performance_solver(T = Float64; grid_size = (48, 48, 24), linear_solver = KrylovJL_GMRES())
    return ConvergentBornSolver(
        wavelength = T(500e-9),           # 500 nm
        permittivity_bg = T(1.33^2),      # Water background
        resolution = (T(50e-9), T(50e-9), T(50e-9)),  # 50 nm isotropic
        grid_size = grid_size,
        boundary_thickness = (T(0.0), T(0.0), T(300e-9)),
        field_attenuation = (T(0.0), T(0.0), T(300e-9)),
        field_attenuation_sharpness = T(1.0),
        periodic_boundary = (true, true, false),
        iterations_max = 10,              # Moderate iterations for performance testing
        tolerance = T(1e-6),              # Realistic tolerance
        linear_solver = linear_solver
    )
end

"""
Create performance test permittivity with multiple scatterers.
"""
function create_performance_permittivity(solver::ConvergentBornSolver{T}) where {T}
    grid_size = solver.grid_size
    eps_bg = solver.permittivity_bg
    permittivity = fill(Complex{T}(eps_bg), grid_size)

    # Add multiple spherical scatterers for more complex scattering
    centers = [
        (grid_size[1]÷3, grid_size[2]÷3, grid_size[3]÷2),
        (2*grid_size[1]÷3, grid_size[2]÷3, grid_size[3]÷2),
        (grid_size[1]÷2, 2*grid_size[2]÷3, grid_size[3]÷3)
    ]

    radius_pixels = 4
    scatterer_permittivity = Complex{T}(1.46^2)  # SiO2

    for center in centers
        for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
            dist_sq = (i - center[1])^2 + (j - center[2])^2 + (k - center[3])^2
            if dist_sq ≤ radius_pixels^2
                permittivity[i, j, k] = scatterer_permittivity
            end
        end
    end

    return permittivity
end

@testset "LinearSolve.jl Performance Tests" begin
    @testset "Solver Algorithm Comparison" begin
        algorithms = [KrylovJL_GMRES(), KrylovJL_BICGSTAB()]

        for algorithm in algorithms
            @testset "Performance: $(typeof(algorithm))" begin
                solver = create_performance_solver(; linear_solver = algorithm)
                permittivity = create_performance_permittivity(solver)

                source = PlaneWaveSource(
                    wavelength = solver.wavelength,
                    polarization = [1.0, 0.0, 0.0],
                    k_vector = [0.0, 0.0, 1.0],
                    amplitude = 1.0
                )

                # Test that solving completes without error
                @test_nowarn E_field, H_field = solve(solver, source, permittivity)

                # Basic timing test (not a strict benchmark, just ensures reasonable performance)
                start_time = time()
                E_field, H_field = solve(solver, source, permittivity)
                elapsed = time() - start_time

                # Should complete in reasonable time (very loose bound for CI)
                @test elapsed < 60.0  # 1 minute max for small problem

                # Test output quality
                @test size(E_field) == (solver.grid_size..., 3)
                @test size(H_field) == (solver.grid_size..., 3)
                @test all(isfinite.(E_field))
                @test all(isfinite.(H_field))
            end
        end
    end

    @testset "Grid Size Scaling" begin
        grid_sizes = [(32, 32, 16), (48, 48, 24)]

        for grid_size in grid_sizes
            @testset "Grid $(grid_size)" begin
                solver = create_performance_solver(; grid_size = grid_size, linear_solver = KrylovJL_GMRES())
                permittivity = create_performance_permittivity(solver)

                source = PlaneWaveSource(
                    wavelength = solver.wavelength,
                    polarization = [1.0, 0.0, 0.0],
                    k_vector = [0.0, 0.0, 1.0],
                    amplitude = 1.0
                )

                # Test that solver works for different grid sizes
                @test_nowarn E_field, H_field = solve(solver, source, permittivity)

                # Verify outputs have correct dimensions
                E_field, H_field = solve(solver, source, permittivity)
                @test size(E_field) == (grid_size..., 3)
                @test size(H_field) == (grid_size..., 3)
            end
        end
    end

    @testset "Memory Usage Validation" begin
        solver = create_performance_solver(; grid_size = (32, 32, 16), linear_solver = KrylovJL_GMRES())
        permittivity = create_performance_permittivity(solver)

        source = PlaneWaveSource(
            wavelength = solver.wavelength,
            polarization = [1.0, 0.0, 0.0],
            k_vector = [0.0, 0.0, 1.0],
            amplitude = 1.0
        )

        # Basic memory usage test - ensure no obvious memory leaks
        gc()  # Force garbage collection
        mem_before = Base.gc_live_bytes()

        for i in 1:3  # Run multiple times
            E_field, H_field = solve(solver, source, permittivity)
        end

        gc()  # Force garbage collection
        mem_after = Base.gc_live_bytes()

        # Memory should not grow excessively (very loose bound)
        mem_growth = mem_after - mem_before
        @test mem_growth < 1_000_000_000  # Less than 1GB growth
    end

    @testset "Precision Validation" begin
        @testset "Float32 Performance" begin
            solver32 = create_performance_solver(
                Float32; grid_size = (24, 24, 12), linear_solver = KrylovJL_GMRES())
            permittivity32 = create_performance_permittivity(solver32)

            source32 = PlaneWaveSource(
                wavelength = solver32.wavelength,
                polarization = [1.0f0, 0.0f0, 0.0f0],
                k_vector = [0.0f0, 0.0f0, 1.0f0],
                amplitude = 1.0f0
            )

            @test_nowarn E_field, H_field = solve(solver32, source32, permittivity32)

            E_field, H_field = solve(solver32, source32, permittivity32)
            @test eltype(E_field) == ComplexF32
            @test eltype(H_field) == ComplexF32
        end
    end
end
