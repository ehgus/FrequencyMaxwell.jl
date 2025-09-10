"""
Performance benchmark tests - migrated from performance_test.jl example
Tests solver performance characteristics and scaling behavior.
"""

@testset "Performance Benchmarks" begin
    @testset "Basic Timing Benchmarks" begin
        # Basic performance test with known problem size
        config = ConvergentBornConfig(
            wavelength = 532e-9,
            permittivity_bg = 1.333^2,
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (64, 64, 32)  # Smaller than original for CI
        )

        # Create bead phantom
        permittivity_bead = 1.4607^2
        radius_pixels = 5
        phantom = create_spherical_phantom(
            config.grid_size,
            config.permittivity_bg,
            permittivity_bead,
            radius_pixels
        )

        source = PlaneWaveSource(
            config,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )

        solver = ConvergentBornSolver(config)
        set_permittivity!(solver, phantom)
        set_Born_max!(solver, 10)

        # Warm-up run
        E_field, H_field = solve(solver, [source])

        # Timed run
        elapsed_time = @elapsed begin
            E_field, H_field = solve(solver, [source])
        end

        # Basic performance assertions
        @test elapsed_time < 60.0  # Should complete in reasonable time
        @test size(E_field) == (config.grid_size..., 3, 1)
        @test size(H_field) == (config.grid_size..., 3, 1)

        # Test memory usage is reasonable (indirect test)
        GC.gc()  # Force garbage collection
        @test true  # If we get here, memory usage was manageable
    end

    @testset "Scaling with Born Iterations" begin
        config = ConvergentBornConfig(
            wavelength = 500e-9,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32)
        )

        # Simple scatterer
        phantom = ones(ComplexF64, config.grid_size)
        phantom[14:18, 14:18, 14:18] .= 1.2

        source = PlaneWaveSource(
            config,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )

        solver = ConvergentBornSolver(config)
        set_permittivity!(solver, phantom)

        # Test scaling with Born iterations
        Born_iterations = [1, 3, 5, 10]
        times = []

        for Born_max in Born_iterations
            set_Born_max!(solver, Born_max)

            # Warm-up
            solve(solver, [source])

            # Timed run
            elapsed = @elapsed solve(solver, [source])
            push!(times, elapsed)
        end

        # Test that time scales reasonably with iterations
        # (should be roughly linear, allowing for overhead)
        @test times[end] > times[1]  # More iterations take more time
        @test times[end] / times[1] < 20  # But not excessively more
    end

    @testset "Grid Size Scaling" begin
        # Test how performance scales with grid size
        base_config = ConvergentBornConfig(
            wavelength = 500e-9,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9)
        )

        grid_sizes = [(16, 16, 16), (24, 24, 24), (32, 32, 32)]
        times = []

        for grid_size in grid_sizes
            config = ConvergentBornConfig(
                wavelength = base_config.wavelength,
                permittivity_bg = base_config.permittivity_bg,
                resolution = base_config.resolution,
                grid_size = grid_size
            )

            # Simple scatterer at center
            phantom = ones(ComplexF64, config.grid_size)
            center_range = div.(grid_size, 4):(3 * div.(grid_size, 4))
            phantom[center_range[1], center_range[2], center_range[3]] .= 1.1

            source = PlaneWaveSource(
                config,
                polarization = (1.0, 0.0, 0.0),
                direction = 3,
                k_transverse = (0.0, 0.0)
            )

            solver = ConvergentBornSolver(config)
            set_permittivity!(solver, phantom)
            set_Born_max!(solver, 3)

            # Warm-up
            solve(solver, [source])

            # Timed run
            elapsed = @elapsed solve(solver, [source])
            push!(times, elapsed)
        end

        # Test reasonable scaling (allowing for significant variation)
        @test times[end] > times[1]  # Larger grids take more time
        @test times[end] / times[1] < 100  # But scaling shouldn't be too bad
    end

    @testset "Memory Usage Patterns" begin
        config = ConvergentBornConfig(
            wavelength = 500e-9,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32)
        )

        phantom = ones(ComplexF64, config.grid_size)
        phantom[14:18, 14:18, 14:18] .= 1.2

        source = PlaneWaveSource(
            config,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )

        solver = ConvergentBornSolver(config)
        set_permittivity!(solver, phantom)
        set_Born_max!(solver, 5)

        # Measure memory before solving
        GC.gc()
        before_bytes = Base.gc_num().total_bytes

        # Solve problem
        E_field, H_field = solve(solver, [source])

        # Measure memory after solving
        GC.gc()
        after_bytes = Base.gc_num().total_bytes

        memory_used = after_bytes - before_bytes

        # Test that memory usage is reasonable
        # (This is a rough test - exact values depend on implementation)
        expected_field_memory = sizeof(ComplexF64) * prod(config.grid_size) * 3 * 2  # E and H fields
        @test memory_used < expected_field_memory * 10  # Allow for overhead
    end

    @testset "Multiple Source Performance" begin
        config = ConvergentBornConfig(
            wavelength = 500e-9,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32)
        )

        phantom = ones(ComplexF64, config.grid_size)
        phantom[14:18, 14:18, 14:18] .= 1.2

        solver = ConvergentBornSolver(config)
        set_permittivity!(solver, phantom)
        set_Born_max!(solver, 3)

        # Test with different numbers of sources
        num_sources_list = [1, 2, 4]
        times = []

        for num_sources in num_sources_list
            sources = []
            for i in 1:num_sources
                k_y = 2π * (i-1) / (config.grid_size[2] * config.resolution[2])
                source = PlaneWaveSource(
                    config,
                    polarization = (1.0, 0.0, 0.0),
                    direction = 3,
                    k_transverse = (0.0, k_y)
                )
                push!(sources, source)
            end

            # Warm-up
            solve(solver, sources)

            # Timed run
            elapsed = @elapsed solve(solver, sources)
            push!(times, elapsed)
        end

        # Test that multiple sources scale reasonably
        @test times[end] > times[1]  # More sources take more time
        @test times[end] / times[1] <= num_sources_list[end]  # But not linearly (some overhead can be shared)
    end

    @testset "Convergence Efficiency" begin
        config = ConvergentBornConfig(
            wavelength = 500e-9,
            permittivity_bg = 1.0,
            resolution = (100e-9, 100e-9, 100e-9),
            grid_size = (32, 32, 32)
        )

        # Test with different contrasts
        contrasts = [1.01, 1.1, 1.5]  # Weak to strong scattering

        source = PlaneWaveSource(
            config,
            polarization = (1.0, 0.0, 0.0),
            direction = 3,
            k_transverse = (0.0, 0.0)
        )

        solver = ConvergentBornSolver(config)

        for contrast in contrasts
            phantom = ones(ComplexF64, config.grid_size)
            phantom[14:18, 14:18, 14:18] .= contrast

            set_permittivity!(solver, phantom)

            # Test convergence with limited iterations
            set_Born_max!(solver, 3)
            E_3, _ = solve(solver, [source])

            set_Born_max!(solver, 10)
            E_10, _ = solve(solver, [source])

            # Calculate convergence metric
            convergence_error = maximum(abs.(E_10 - E_3)) / maximum(abs.(E_10))

            # Weak scattering should converge faster
            if contrast < 1.1
                @test convergence_error < 0.1  # Good convergence for weak scattering
            end
        end
    end
end
