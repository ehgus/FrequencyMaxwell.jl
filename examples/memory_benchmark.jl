"""
Memory Benchmarking Suite for FrequencyMaxwell Multi-Source Electromagnetic Solver

This comprehensive benchmark measures memory usage patterns, scaling behavior, and
performance bottlenecks in FrequencyMaxwell's multi-source electromagnetic solver
to guide memory optimization efforts.

# Key Features:
- Memory scaling analysis with source count (1, 2, 4, 8, 16+ sources)
- Grid size scaling benchmarks (64³, 128³, 256³)
- Peak memory usage tracking and allocation patterns
- Garbage collection impact measurement
- Memory stress testing to identify bottlenecks
- CSV output for optimization tracking
- Visualization of memory scaling patterns

# Usage:
```julia
julia examples/memory_benchmark.jl
```

This benchmark follows the basic_scattering.jl structure and integrates with
the refactored multi-source API for realistic memory profiling.
"""

using FrequencyMaxwell
using LinearSolve
using LinearAlgebra
using Printf
using Statistics
using Dates

# Try to use CSV and DataFrames for data output, fallback gracefully  
CSV_AVAILABLE = try
    using CSV
    using DataFrames
    true
catch
    @warn "CSV/DataFrames not available, skipping CSV output"
    false
end

# Try to use WGLMakie for visualization, fallback gracefully
PLOTTING_AVAILABLE = try
    using WGLMakie
    true
catch
    @warn "WGLMakie not available, skipping visualization components"
    false
end

"""
    MemoryStats
    
Structure to hold comprehensive memory measurement results.
"""
struct MemoryStats
    # Basic measurements
    peak_memory_mb::Float64
    total_allocations::Int64
    allocation_bytes::Int64
    gc_time_ms::Float64
    solve_time_s::Float64

    # Scaling metrics
    sources_count::Int
    grid_size::Tuple{Int, Int, Int}
    memory_per_source_mb::Float64
    memory_per_voxel_bytes::Float64

    # Performance metrics
    memory_efficiency::Float64  # Theoretical vs actual memory usage
    gc_pressure::Float64       # GC time / total time
end

"""
    BenchmarkConfig
    
Configuration for memory benchmarking scenarios.
"""
struct BenchmarkConfig
    # Base electromagnetic parameters (following basic_scattering.jl)
    wavelength::Float64
    permittivity_bg::Float64
    resolution::Tuple{Float64, Float64, Float64}

    # Solver parameters
    tolerance::Float64
    iterations_max::Int

    # Benchmark parameters
    max_sources::Int
    grid_sizes::Vector{Tuple{Int, Int, Int}}
    stress_test_enabled::Bool
    warmup_runs::Int
    measurement_runs::Int

    function BenchmarkConfig(;
            wavelength = 532e-9,                    # 532 nm green laser
            permittivity_bg = 1.333^2,             # Water background
            resolution = (50e-9, 50e-9, 50e-9),   # 50 nm isotropic
            tolerance = 1e-6,
            iterations_max = -1,                   # Auto-calculate
            max_sources = 16,
            grid_sizes = [(64, 64, 32), (128, 128, 64), (256, 256, 128)],
            stress_test_enabled = true,
            warmup_runs = 2,
            measurement_runs = 3
    )
        new(wavelength, permittivity_bg, resolution, tolerance, iterations_max,
            max_sources, grid_sizes, stress_test_enabled, warmup_runs, measurement_runs)
    end
end

"""
    create_solver_config(benchmark_config, grid_size)
    
Create ConvergentBornConfig from benchmark parameters.
"""
function create_solver_config(benchmark_config::BenchmarkConfig, grid_size::Tuple{
        Int, Int, Int})
    return ConvergentBornConfig(
        wavelength = benchmark_config.wavelength,
        permittivity_bg = benchmark_config.permittivity_bg,
        resolution = benchmark_config.resolution,
        grid_size = grid_size,
        boundary_thickness = (0.0, 0.0, 3.0e-6),      # 3 μm padding in Z
        field_attenuation = (0.0, 0.0, 3.0e-6),       # 3 μm attenuation in Z
        field_attenuation_sharpness = 1.0,
        periodic_boundary = (true, true, false),        # Periodic in XY, absorbing in Z
        iterations_max = benchmark_config.iterations_max,
        tolerance = benchmark_config.tolerance,
        linear_solver = KrylovJL_GMRES()               # Direct LinearSolver object
    )
end

"""
    create_test_phantom(grid_size)
    
Create a realistic test phantom (spherical bead) for benchmarking.
"""
function create_test_phantom(grid_size::Tuple{Int, Int, Int})
    # Create spherical SiO2 bead (similar to basic_scattering.jl)
    bead_radius_pixels = min(10.0, minimum(grid_size) / 8)  # Scale with grid size

    return phantom_bead(
        grid_size,
        [1.46^2],           # SiO2 bead (n=1.46)
        bead_radius_pixels
    )
end

"""
    create_multi_source_array(n_sources, wavelength)
    
Create array of coherent plane wave sources for multi-source testing.
"""
function create_multi_source_array(n_sources::Int, wavelength::Float64)
    sources = PlaneWaveSource{Float64}[]

    if n_sources == 1
        # Single source baseline
        push!(sources,
            PlaneWaveSource(
                wavelength = wavelength,
                polarization = [1.0, 0.0, 0.0],    # X-polarized
                k_vector = [0.0, 0.0, 1.0],        # +Z propagation
                amplitude = 1.0
            ))
    else
        # Multi-source configuration: create interference pattern
        # Use angular distribution to create realistic multi-source scenario
        for i in 1:n_sources
            # Create sources with slight angular variations for interference
            angle = 2π * (i - 1) / n_sources  # Distribute around circle

            # Small angular deviation (±5 degrees max)
            θ = 0.1 * sin(angle)  # ~5.7 degrees max deviation
            φ = 0.1 * cos(angle)

            k_x = sin(θ) * cos(φ)
            k_y = sin(θ) * sin(φ)
            k_z = cos(θ)

            # Normalize k-vector
            k_vec = [k_x, k_y, k_z]
            k_hat = k_vec ./ norm(k_vec)

            # Create polarization perpendicular to k-vector
            # Start with [1,0,0] and make it perpendicular
            pol_base = [1.0, 0.0, 0.0]
            pol_perp = pol_base .- (dot(pol_base, k_hat) .* k_hat)
            pol_perp = pol_perp ./ norm(pol_perp)  # Normalize

            push!(sources,
                PlaneWaveSource(
                    wavelength = wavelength,
                    polarization = pol_perp,                   # Transverse polarization
                    k_vector = k_hat,                          # Unit k-vector
                    amplitude = 1.0 / sqrt(n_sources),        # Normalize total power
                    phase = 2π * rand()                        # Random phase for realism
                ))
        end
    end

    return sources
end

"""
    measure_memory_usage(func)
    
Measure comprehensive memory usage statistics for a function call.
"""
function measure_memory_usage(func::Function)
    # Force garbage collection before measurement
    GC.gc()

    # Get initial memory state
    initial_memory = Sys.total_memory() - Sys.free_memory()

    # Measure execution with detailed allocation tracking
    stats = @timed func()
    solve_time = stats.time
    allocation_bytes = stats.bytes
    gc_time = stats.gctime

    # Force GC again and measure peak
    GC.gc()

    # Get final memory state
    final_memory = Sys.total_memory() - Sys.free_memory()
    peak_memory = max(final_memory, initial_memory + allocation_bytes)

    return (
        peak_memory_mb = peak_memory / (1024^2),
        allocation_bytes = allocation_bytes,
        gc_time_ms = gc_time * 1000,
        solve_time_s = solve_time
    ),
    stats.value
end

"""
    benchmark_single_scenario(config, grid_size, n_sources)
    
Benchmark a single scenario with specified grid size and source count.
"""
function benchmark_single_scenario(
        benchmark_config::BenchmarkConfig,
        grid_size::Tuple{Int, Int, Int},
        n_sources::Int
)
    @printf("  Grid: %dx%dx%d, Sources: %d\n", grid_size..., n_sources)

    # Create solver configuration
    solver_config = create_solver_config(benchmark_config, grid_size)

    # Create phantom
    phantom = create_test_phantom(grid_size)

    # Create sources
    sources = create_multi_source_array(n_sources, benchmark_config.wavelength)

    # Warmup runs (not measured, for JIT compilation)
    @printf("    Warming up...")
    for _ in 1:benchmark_config.warmup_runs
        try
            solver = ConvergentBornSolver(solver_config)
            E_field, H_field = solve(solver, sources, phantom)
        catch e
            @warn "Warmup failed: $e"
            break
        end
    end
    @printf(" done\n")

    # Measurement runs
    memory_measurements = []
    solve_times = []

    @printf("    Measuring...")
    for run in 1:benchmark_config.measurement_runs
        try
            memory_stats, result = measure_memory_usage() do
                solver = ConvergentBornSolver(solver_config)
                return solve(solver, sources, phantom)
            end

            push!(memory_measurements, memory_stats)
            push!(solve_times, memory_stats.solve_time_s)

            # Validate result
            E_field, H_field = result
            if any(isnan.(E_field)) || any(isnan.(H_field))
                @warn "NaN detected in solution for run $run"
            end

        catch e
            @error "Measurement run $run failed: $e"
            # Return early if measurements are failing
            if length(memory_measurements) == 0
                return nothing
            end
            break
        end
    end
    @printf(" done\n")

    # Compute statistics from multiple runs
    if isempty(memory_measurements)
        @warn "No successful measurements for grid $grid_size, sources $n_sources"
        return nothing
    end

    # Average the measurements
    avg_peak_memory = mean(m.peak_memory_mb for m in memory_measurements)
    avg_allocation_bytes = mean(m.allocation_bytes for m in memory_measurements)
    avg_gc_time = mean(m.gc_time_ms for m in memory_measurements)
    avg_solve_time = mean(m.solve_time_s for m in memory_measurements)

    # Calculate derived metrics
    total_voxels = prod(grid_size)
    memory_per_source = avg_peak_memory / n_sources
    memory_per_voxel = avg_allocation_bytes / total_voxels
    gc_pressure = avg_gc_time / (avg_solve_time * 1000)  # GC fraction of total time

    # Estimate theoretical memory usage (rough calculation)
    # Complex{Float64} E and H fields = 16 bytes per component * 3 components * 2 fields
    theoretical_field_memory = total_voxels * 16 * 3 * 2 / (1024^2)  # MB
    memory_efficiency = theoretical_field_memory / avg_peak_memory

    return MemoryStats(
        avg_peak_memory,
        Int64(round(avg_allocation_bytes / 16)),  # Rough allocation count
        Int64(round(avg_allocation_bytes)),
        avg_gc_time,
        avg_solve_time,
        n_sources,
        grid_size,
        memory_per_source,
        memory_per_voxel,
        memory_efficiency,
        gc_pressure
    )
end

"""
    run_source_scaling_benchmark(config)
    
Test memory scaling with increasing number of sources.
"""
function run_source_scaling_benchmark(benchmark_config::BenchmarkConfig)
    println("\n" * "="^60)
    println("SOURCE SCALING BENCHMARK")
    println("="^60)

    # Use medium grid size for source scaling tests
    grid_size = benchmark_config.grid_sizes[2]  # Typically (128,128,64)

    source_counts = [1, 2, 4, 8]
    if benchmark_config.max_sources >= 16
        push!(source_counts, 16)
    end

    results = MemoryStats[]

    for n_sources in source_counts
        @printf("\nBenchmarking %d source(s):\n", n_sources)

        result = benchmark_single_scenario(benchmark_config, grid_size, n_sources)
        if result !== nothing
            push!(results, result)

            # Print summary for this test
            @printf("  Peak Memory: %.1f MB (%.2f MB/source)\n",
                result.peak_memory_mb, result.memory_per_source_mb)
            @printf("  Solve Time: %.2f s\n", result.solve_time_s)
            @printf("  GC Pressure: %.1f%%\n", result.gc_pressure * 100)
        end
    end

    return results
end

"""
    run_grid_scaling_benchmark(config)
    
Test memory scaling with increasing grid sizes.
"""
function run_grid_scaling_benchmark(benchmark_config::BenchmarkConfig)
    println("\n" * "="^60)
    println("GRID SIZE SCALING BENCHMARK")
    println("="^60)

    results = MemoryStats[]

    for grid_size in benchmark_config.grid_sizes
        @printf("\nBenchmarking grid %dx%dx%d:\n", grid_size...)

        # Test with 4 sources for realistic multi-source scenario
        n_sources = min(4, benchmark_config.max_sources)

        result = benchmark_single_scenario(benchmark_config, grid_size, n_sources)
        if result !== nothing
            push!(results, result)

            # Print summary
            total_voxels = prod(grid_size)
            @printf("  Peak Memory: %.1f MB (%.2f bytes/voxel)\n",
                result.peak_memory_mb, result.memory_per_voxel_bytes)
            @printf("  Solve Time: %.2f s\n", result.solve_time_s)
            @printf("  Memory Efficiency: %.1f%%\n", result.memory_efficiency * 100)
        end
    end

    return results
end

"""
    run_stress_test(config)
    
Run stress tests to identify memory limits and bottlenecks.
"""
function run_stress_test(benchmark_config::BenchmarkConfig)
    if !benchmark_config.stress_test_enabled
        return MemoryStats[]
    end

    println("\n" * "="^60)
    println("MEMORY STRESS TEST")
    println("="^60)

    results = MemoryStats[]

    # Test with largest grid and maximum sources
    largest_grid = benchmark_config.grid_sizes[end]
    max_sources = benchmark_config.max_sources

    @printf("\nStress testing with grid %dx%dx%d and %d sources:\n",
        largest_grid..., max_sources)

    result = benchmark_single_scenario(benchmark_config, largest_grid, max_sources)
    if result !== nothing
        push!(results, result)

        @printf("  STRESS TEST RESULTS:\n")
        @printf("  Peak Memory: %.1f MB\n", result.peak_memory_mb)
        @printf("  Memory per Source: %.2f MB\n", result.memory_per_source_mb)
        @printf("  Solve Time: %.2f s\n", result.solve_time_s)
        @printf("  GC Pressure: %.1f%%\n", result.gc_pressure * 100)

        # Memory stress analysis
        if result.peak_memory_mb > 8000  # > 8 GB
            @warn "High memory usage detected ($(result.peak_memory_mb) MB)"
        end

        if result.gc_pressure > 0.2  # > 20% time in GC
            @warn "High GC pressure detected ($(result.gc_pressure*100)%)"
        end

        if result.memory_per_source_mb > 500  # > 500 MB per source
            @warn "Memory scaling may be suboptimal ($(result.memory_per_source_mb) MB/source)"
        end

    else
        @error "Stress test failed - unable to complete benchmark"
    end

    return results
end

"""
    save_results_csv(results, filename)
    
Save benchmark results to CSV file for tracking optimization progress.
"""
function save_results_csv(results::Vector{Vector{MemoryStats}}, filename::String)
    if !CSV_AVAILABLE
        @warn "CSV functionality not available, skipping file output"
        return
    end

    if isempty(results) || all(isempty, results)
        @warn "No results to save"
        return
    end

    # Flatten results
    all_results = vcat(results...)

    # Create DataFrame
    df = DataFrame(
        timestamp = fill(now(), length(all_results)),
        test_type = String[],
        sources_count = [r.sources_count for r in all_results],
        grid_x = [r.grid_size[1] for r in all_results],
        grid_y = [r.grid_size[2] for r in all_results],
        grid_z = [r.grid_size[3] for r in all_results],
        peak_memory_mb = [r.peak_memory_mb for r in all_results],
        allocation_bytes = [r.allocation_bytes for r in all_results],
        gc_time_ms = [r.gc_time_ms for r in all_results],
        solve_time_s = [r.solve_time_s for r in all_results],
        memory_per_source_mb = [r.memory_per_source_mb for r in all_results],
        memory_per_voxel_bytes = [r.memory_per_voxel_bytes for r in all_results],
        memory_efficiency = [r.memory_efficiency for r in all_results],
        gc_pressure = [r.gc_pressure for r in all_results]
    )

    # Add test type classification
    test_types = String[]
    result_idx = 1
    for (i, result_group) in enumerate(results)
        if i == 1
            append!(test_types, fill("source_scaling", length(result_group)))
        elseif i == 2
            append!(test_types, fill("grid_scaling", length(result_group)))
        else
            append!(test_types, fill("stress_test", length(result_group)))
        end
    end
    df.test_type = test_types

    # Save to CSV
    CSV.write(filename, df)
    @printf("Results saved to %s\n", filename)
end

"""
    create_memory_plots(results)
    
Create visualization plots of memory scaling patterns.
"""
function create_memory_plots(results::Vector{Vector{MemoryStats}})
    if !PLOTTING_AVAILABLE
        @info "Plotting not available, skipping visualization"
        return
    end

    if isempty(results) || all(isempty, results)
        @warn "No results to plot"
        return
    end

    try
        # Source scaling plot
        if !isempty(results[1])
            source_results = results[1]
            sources = [r.sources_count for r in source_results]
            memory = [r.peak_memory_mb for r in source_results]

            fig1 = Figure(resolution = (800, 600))
            ax1 = Axis(fig1[1, 1],
                xlabel = "Number of Sources",
                ylabel = "Peak Memory (MB)",
                title = "Memory Scaling vs Source Count")

            lines!(ax1, sources, memory, marker = :circle, linewidth = 3)

            # Add ideal linear scaling reference
            if length(sources) > 1
                linear_ref = memory[1] .* (sources ./ sources[1])
                lines!(ax1, sources, linear_ref, linestyle = :dash,
                    color = :red, label = "Linear scaling")
                axislegend(ax1)
            end

            display(fig1)
        end

        # Grid scaling plot
        if length(results) > 1 && !isempty(results[2])
            grid_results = results[2]
            voxels = [prod(r.grid_size) for r in grid_results]
            memory = [r.peak_memory_mb for r in grid_results]

            fig2 = Figure(resolution = (800, 600))
            ax2 = Axis(fig2[1, 1],
                xlabel = "Total Voxels",
                ylabel = "Peak Memory (MB)",
                title = "Memory Scaling vs Grid Size",
                xscale = log10, yscale = log10)

            scatter!(ax2, voxels, memory, markersize = 15)

            # Fit and display scaling trend
            if length(voxels) > 1
                log_voxels = log10.(voxels)
                log_memory = log10.(memory)

                # Simple linear fit in log-log space
                A = [ones(length(log_voxels)) log_voxels]
                coeffs = A \ log_memory
                scaling_exponent = coeffs[2]

                # Plot fit line
                fit_x = 10 .^ range(minimum(log_voxels), maximum(log_voxels), 100)
                fit_y = 10 .^ (coeffs[1] .+ coeffs[2] .* log10.(fit_x))
                lines!(ax2, fit_x, fit_y, color = :red,
                    label = @sprintf("Scaling: N^%.2f", scaling_exponent))

                axislegend(ax2)

                @printf("Memory scales as N^%.2f with grid size\n", scaling_exponent)
            end

            display(fig2)
        end

    catch e
        @warn "Plotting failed: $e"
    end
end

"""
    print_summary_report(results)
    
Print comprehensive summary of benchmark results.
"""
function print_summary_report(results::Vector{Vector{MemoryStats}})
    println("\n" * "="^60)
    println("MEMORY BENCHMARK SUMMARY REPORT")
    println("="^60)

    if isempty(results) || all(isempty, results)
        println("No benchmark results available.")
        return
    end

    all_results = vcat(results...)

    if isempty(all_results)
        println("No successful benchmark runs completed.")
        return
    end

    # Overall statistics
    println("\nOVERALL STATISTICS:")
    @printf("  Tests completed: %d\n", length(all_results))
    @printf("  Peak memory range: %.1f - %.1f MB\n",
        minimum(r.peak_memory_mb for r in all_results),
        maximum(r.peak_memory_mb for r in all_results))
    @printf("  Solve time range: %.2f - %.2f s\n",
        minimum(r.solve_time_s for r in all_results),
        maximum(r.solve_time_s for r in all_results))

    # Memory efficiency analysis
    avg_efficiency = mean(r.memory_efficiency for r in all_results)
    @printf("  Average memory efficiency: %.1f%%\n", avg_efficiency * 100)

    if avg_efficiency < 0.1  # < 10% efficient
        println("  ⚠️  LOW MEMORY EFFICIENCY - significant optimization potential")
    elseif avg_efficiency < 0.3  # < 30% efficient
        println("  ⚠️  MODERATE MEMORY EFFICIENCY - optimization recommended")
    else
        println("  ✓  REASONABLE MEMORY EFFICIENCY")
    end

    # GC pressure analysis
    avg_gc_pressure = mean(r.gc_pressure for r in all_results)
    @printf("  Average GC pressure: %.1f%%\n", avg_gc_pressure * 100)

    if avg_gc_pressure > 0.2  # > 20%
        println("  ⚠️  HIGH GC PRESSURE - memory allocation optimization needed")
    else
        println("  ✓  REASONABLE GC PRESSURE")
    end

    # Scaling analysis
    println("\nSCALING ANALYSIS:")

    if !isempty(results[1])  # Source scaling results
        source_results = results[1]
        if length(source_results) > 1
            memory_1_source = source_results[1].peak_memory_mb
            memory_max_source = source_results[end].peak_memory_mb
            sources_max = source_results[end].sources_count

            scaling_factor = memory_max_source / memory_1_source
            ideal_scaling = sources_max
            efficiency_ratio = ideal_scaling / scaling_factor

            @printf("  Source scaling: %dx memory for %dx sources (efficiency: %.1f%%)\n",
                round(Int, scaling_factor), sources_max, efficiency_ratio * 100)

            if efficiency_ratio > 0.8  # > 80% efficient scaling
                println("  ✓  EXCELLENT source scaling efficiency")
            elseif efficiency_ratio > 0.5  # > 50% efficient scaling
                println("  ⚠️  MODERATE source scaling efficiency")
            else
                println("  ❌  POOR source scaling efficiency - major optimization needed")
            end
        end
    end

    # Recommendations
    println("\nOPTIMIZATION RECOMMENDATIONS:")

    worst_gc = maximum(r.gc_pressure for r in all_results)
    if worst_gc > 0.3
        println("  1. CRITICAL: Reduce memory allocations (GC pressure: $(round(worst_gc*100))%)")
    end

    if avg_efficiency < 0.2
        println("  2. HIGH: Improve memory utilization efficiency")
    end

    largest_memory = maximum(r.peak_memory_mb for r in all_results)
    if largest_memory > 4000  # > 4 GB
        println("  3. MEDIUM: Consider memory streaming for large problems")
    end

    # Check for sublinear scaling
    if !isempty(results[1]) && length(results[1]) > 1
        source_results = results[1]
        memory_per_source = [r.memory_per_source_mb for r in source_results]
        if maximum(memory_per_source) / minimum(memory_per_source) > 1.5
            println("  4. MEDIUM: Investigate source scaling bottlenecks")
        end
    end

    println("\n" * "="^60)
    println("For 10x larger problems, focus on top recommendations first.")
    println("="^60)
end

"""
    main()
    
Main benchmark execution function.
"""
function main()
    println("FrequencyMaxwell Memory Benchmark Suite")
    println("=" ^ 45)
    println("Target: Memory optimization for 10x larger electromagnetic problems")
    println()

    # Configuration
    benchmark_config = BenchmarkConfig()

    println("Benchmark Configuration:")
    @printf("  Wavelength: %.0f nm\n", benchmark_config.wavelength * 1e9)
    @printf("  Background: n=%.3f\n", sqrt(benchmark_config.permittivity_bg))
    @printf("  Resolution: %.0f nm isotropic\n", benchmark_config.resolution[1] * 1e9)
    @printf("  Max sources: %d\n", benchmark_config.max_sources)
    @printf("  Grid sizes: %s\n",
        join(["$(g[1])×$(g[2])×$(g[3])" for g in benchmark_config.grid_sizes], ", "))
    @printf("  Measurement runs: %d\n", benchmark_config.measurement_runs)

    # Run benchmark suite
    results = Vector{MemoryStats}[]

    try
        # 1. Source scaling benchmark
        source_results = run_source_scaling_benchmark(benchmark_config)
        push!(results, source_results)

        # 2. Grid scaling benchmark  
        grid_results = run_grid_scaling_benchmark(benchmark_config)
        push!(results, grid_results)

        # 3. Stress test
        stress_results = run_stress_test(benchmark_config)
        push!(results, stress_results)

        # Save results
        timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
        csv_filename = "memory_benchmark_results_$(timestamp).csv"
        save_results_csv(results, csv_filename)

        # Create visualizations
        create_memory_plots(results)

        # Print comprehensive report
        print_summary_report(results)

        println("\n✅ Memory benchmark completed successfully!")
        println("📊 Results saved to: $csv_filename")

    catch e
        @error "Benchmark failed with error: $e"
        println("\nStack trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        return false
    end

    return true
end

# Run benchmark if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    success = main()
    exit(success ? 0 : 1)
end
