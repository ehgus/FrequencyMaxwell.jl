# Memory Benchmark Suite for FrequencyMaxwell

## Overview

The `memory_benchmark.jl` script provides comprehensive memory profiling and performance analysis for FrequencyMaxwell's multi-source electromagnetic solver. This benchmark is designed to guide memory optimization efforts and help achieve the goal of handling 10x larger electromagnetic problems.

## Key Features

### Memory Metrics Measured
- **Peak Memory Usage**: Maximum RSS during solve operations
- **Memory Scaling**: How memory grows with source count (1, 2, 4, 8, 16+ sources)
- **Allocation Patterns**: Track allocations during different solve phases
- **Garbage Collection Impact**: GC pressure and timing
- **Memory Efficiency**: Memory per source (should be sublinear)
- **Memory per Voxel**: Scaling with problem size

### Benchmark Scenarios
- **Single Source Baseline**: Current basic_scattering.jl behavior for comparison
- **Multi-Source Scaling**: Test 2, 4, 8, 16 coherent sources with same phantom
- **Grid Size Scaling**: Test different grid sizes (64³, 128³, 256³) to find bottlenecks
- **Memory Stress Test**: Push towards memory limits to identify breaking points

## Usage

### Basic Usage

Run the complete benchmark suite:
```bash
julia --project=. examples/memory_benchmark.jl
```

### Programmatic Usage

```julia
include("examples/memory_benchmark.jl")

# Create custom configuration
config = BenchmarkConfig(
    max_sources = 8,
    grid_sizes = [(64,64,32), (128,128,64)],
    measurement_runs = 3,
    warmup_runs = 2,
    stress_test_enabled = true
)

# Run specific benchmark
results = run_source_scaling_benchmark(config)
```

### Configuration Options

```julia
BenchmarkConfig(;
    wavelength = 532e-9,                    # 532 nm green laser
    permittivity_bg = 1.333^2,             # Water background
    resolution = (50e-9, 50e-9, 50e-9),   # 50 nm isotropic
    tolerance = 1e-6,
    iterations_max = -1,                   # Auto-calculate
    max_sources = 16,
    grid_sizes = [(64,64,32), (128,128,64), (256,256,128)],
    stress_test_enabled = true,
    warmup_runs = 2,
    measurement_runs = 3
)
```

## Output and Analysis

### Console Output
The benchmark provides detailed console output including:
- Real-time progress updates
- Memory usage statistics
- Performance metrics
- Scaling analysis
- Optimization recommendations

### CSV Output (Optional)
If CSV.jl and DataFrames.jl are available, results are saved to timestamped CSV files:
```
memory_benchmark_results_20231204_143022.csv
```

### Visualization (Optional)
If WGLMakie.jl is available, creates plots showing:
- Memory scaling vs source count
- Memory scaling vs grid size
- Scaling exponent analysis

## Understanding Results

### Memory Efficiency
- **> 30%**: Reasonable efficiency
- **10-30%**: Moderate efficiency, optimization recommended
- **< 10%**: Low efficiency, significant optimization potential

### GC Pressure
- **< 10%**: Excellent GC behavior
- **10-20%**: Reasonable GC pressure
- **> 20%**: High GC pressure, allocation optimization needed

### Scaling Analysis
- **Source scaling**: How memory grows with number of sources
  - Ideal: Linear scaling (2x sources = 2x memory)
  - Good: < 1.5x overhead per doubling
  - Poor: > 2x overhead per doubling

- **Grid scaling**: How memory grows with problem size
  - Expected: ~O(N) for N voxels (field storage)
  - Actual: May show higher scaling due to algorithm overhead

## Interpreting Recommendations

### Critical Issues (Address First)
- **High GC Pressure (>30%)**: Indicates excessive allocations
  - Focus on reducing temporary array creation
  - Implement in-place operations
  - Pre-allocate working arrays

### High Priority Issues
- **Low Memory Efficiency (<20%)**: Memory usage far exceeds theoretical minimum
  - Review data structures for redundancy
  - Optimize Green's function storage
  - Consider memory streaming for large arrays

### Medium Priority Issues
- **Large Memory Usage (>4GB)**: Consider memory streaming
- **Poor Source Scaling**: Investigate multi-source bottlenecks

## Integration with Development Workflow

### Before Optimization
1. Run benchmark to establish baseline
2. Identify worst bottlenecks from recommendations
3. Save CSV results for comparison

### After Optimization
1. Run same benchmark configuration
2. Compare CSV results to quantify improvement
3. Verify optimization didn't break functionality

### Target Metrics for 10x Larger Problems
- Memory efficiency: > 20%
- GC pressure: < 10%
- Source scaling efficiency: > 80%
- Peak memory: manageable within system limits

## Technical Details

### Memory Measurement
- Uses Julia's `@timed` macro for accurate allocation tracking
- Forces garbage collection before/after measurements
- Accounts for JIT compilation with warmup runs
- Multiple measurement runs for statistical reliability

### Test Scenarios
- **Phantom**: Spherical SiO2 bead (n=1.46) scaled to grid size
- **Sources**: Coherent plane waves with proper transverse polarization
- **Solver**: ConvergentBornSolver with realistic CBS parameters
- **Physics**: 532nm wavelength in water background (n=1.333)

### Error Handling
- Graceful fallback when optional dependencies unavailable
- Robust error handling for failed solve operations
- Validation of electromagnetic field solutions (NaN detection)

## Dependencies

### Required
- FrequencyMaxwell.jl
- LinearAlgebra (stdlib)
- Printf (stdlib)
- Statistics (stdlib)
- Dates (stdlib)

### Optional (with graceful fallback)
- CSV.jl + DataFrames.jl (for CSV output)
- WGLMakie.jl (for visualization)

## Examples

### Quick Memory Check
```julia
# Test current memory usage with 4 sources on medium grid
result = benchmark_single_scenario(
    BenchmarkConfig(measurement_runs=1),
    (128,128,64),
    4
)
println("Memory per source: $(result.memory_per_source_mb) MB")
```

### Scaling Analysis
```julia
# Compare memory scaling between different source counts
for n_sources in [1, 2, 4, 8]
    result = benchmark_single_scenario(config, (128,128,64), n_sources)
    efficiency = n_sources / (result.peak_memory_mb / result.memory_per_source_mb)
    println("$n_sources sources: $(round(efficiency*100))% scaling efficiency")
end
```

This benchmark provides the foundation for systematic memory optimization of FrequencyMaxwell's multi-source electromagnetic solver.