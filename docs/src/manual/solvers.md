# Solvers

FrequencyMaxwell.jl implements sophisticated electromagnetic solvers for frequency-domain Maxwell equations. The current implementation focuses on the Convergent Born Iterative Method (CBS), which has been reformulated to leverage modern linear algebra frameworks.

## Convergent Born Iterative Solver

The `ConvergentBornSolver` implements the Convergent Born Iterative Method, a powerful technique for solving electromagnetic scattering problems in complex media.

### Mathematical Foundation

The CBS method solves the frequency-domain Maxwell equations by reformulating the scattering problem as a linear system:

```math
(\mathbf{I} - \mathbf{G} \cdot \mathbf{V}) \mathbf{E} = \mathbf{G} \cdot \mathbf{E}_{\text{inc}}
```

Where:
- ``\mathbf{I}`` is the identity operator
- ``\mathbf{G}`` is the dyadic Green's function operator
- ``\mathbf{V}`` represents the permittivity contrast
- ``\mathbf{J}`` is the induced current density
- ``\mathbf{E}^{\text{inc}}`` is the incident electric field

### Solver Architecture

#### ConvergentBornSolver Type

```julia
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32)
)
```

The solver type integrates:
- **Configuration fields**: Physical and numerical parameters
- **Mutable state**: Arrays for fields, Green's functions, and working memory
- **LinearSolve.jl integration**: Custom operators for flexible solving

#### Integrated Design

The solver combines configuration and state in a single object:

- **Direct field access**: `solver.wavelength`, `solver.grid_size`, etc.
- **Type-safe validation**: All parameters validated during construction
- **Utility functions**: `domain_size(solver)`, `wavenumber_background(solver)`

This design provides:
- Simplified API with single solver object
- Direct access to all configuration parameters

### Core Solving Function

The main interface is the `solve` function, which uses multiple dispatch to handle different source types:

```julia
# Single source
Efield, Hfield = solve(solver, source, permittivity_phantom)

# Multiple sources (coherent superposition)
Efield, Hfield = solve(solver, [source1, source2], permittivity_phantom)
```

### LinearSolve.jl Integration

FrequencyMaxwell.jl leverages the LinearSolve.jl ecosystem for maximum flexibility and performance.

#### Custom Linear Operators

**CBSLinearOperator**

Implements the $(I - G \cdot V)$ operator as a LinearSolve.jl-compatible linear operator:

```julia
# The operator is applied as: (I - G*V) * J
# Where J is the current density vector
```

This reformulation enables:
- Use of any LinearSolve.jl algorithm (GMRES, BiCGStab, etc.)
- Automatic algorithm selection based on problem characteristics
- Integration with the broader SciML ecosystem


#### Supported Algorithms

Any LinearSolve.jl algorithm can be used. Tested algorithms include:

- **KrylovJL_GMRES()**: Recommended default Krylov method
  - Optimized implementation from Krylov.jl
  - Good general-purpose choice with excellent performance
  - Memory efficient with restart capabilities

- **GMRES()**: Standard Generalized Minimal RESidual method
  - Alternative implementation
  - Memory requirements scale with restart parameter

- **BiCGStabL(l)**: Bi-Conjugate Gradient Stabilized
  - Often faster convergence than GMRES
  - Parameter `l` controls stability (typically 2-4)

- **IDRsAlgorithm(s)**: Induced Dimension Reduction
  - Can be very fast for certain problems
  - Parameter `s` controls memory vs. speed tradeoff

Example configuration:

```julia
using LinearSolve

config = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64),
    linear_solver = KrylovJL_GMRES(),  # Recommended default algorithm
    tolerance = 1e-6
)
```

### Multi-Source Support

The solver natively supports multiple coherent electromagnetic sources:

```julia
# Define multiple sources
source1 = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [1.0, 0.0, 0.0],
    k_vector = [0.0, 0.0, 1.0],
    phase = 0.0
)

source2 = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [0.0, 1.0, 0.0],
    k_vector = [0.0, 0.0, 1.0],
    phase = π/2  # 90° phase shift
)

# Solve with coherent superposition
Efield, Hfield = solve(solver, [source1, source2], phantom)
```

This enables simulation of:
- Interference patterns
- Polarization control
- Complex illumination schemes
- Multi-beam configurations

### Precision Selection

Choose precision based on problem requirements:

```julia
# High precision (Float64) - default
config64 = ConvergentBornConfig{Float64}(...)

# Memory-efficient (Float32)
config32 = ConvergentBornConfig{Float32}(...)
```

Float32 provides:
- 50% memory reduction
- Faster computation on many GPUs
- Usually sufficient precision for most applications

### Memory Management

#### Preallocated Arrays

The solver preallocates all major arrays during initialization:
- Field arrays for E and H fields
- Working arrays for iterative computation
- FFT planning for optimized Green's function evaluation

The solver uses in-place operations to minimize allocations and includes garbage collection hints for large arrays.

### Performance Characteristics

#### Computational Complexity

The CBS solver has complexity characteristics:
- **Memory**: O(N³) for N³ grid points
- **Time per iteration**: O(N³ log N) due to FFT operations
- **Iterations to convergence**: Problem-dependent, typically 10-100

#### Convergence Behavior

Convergence depends on:
- **Permittivity contrast**: Higher contrast requires more iterations
- **Grid resolution**: Finer grids may require more iterations
- **Linear solver choice**: Different algorithms have different convergence rates

### Convergence Monitoring and Diagnostics

#### Automatic Convergence Detection

The solver automatically monitors convergence using the LinearSolve.jl framework. Convergence information is available through the linear solver's standard interface.

#### Convergence Troubleshooting

If convergence issues occur:

1. **Try different linear solvers**:
```julia
# Try BiCGStabL instead of KrylovJL_GMRES
config = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64),
    linear_solver = BiCGStabL(2)  # Alternative algorithm
)
```

2. **Adjust tolerance**:
```julia
config = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64),
    tolerance = 1e-5  # Relax from default 1e-6
)
```

3. **Try different grid resolution**:
```julia
config = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (75e-9, 75e-9, 75e-9),  # Coarser resolution
    grid_size = (96, 96, 48),             # Smaller problem size
    linear_solver = KrylovJL_GMRES()
)
```

### Advanced Features

#### State Reuse

Solver state can be reused across multiple solves:

```julia
solver = ConvergentBornSolver(config)

# Multiple solves reuse allocated memory
for phantom in phantom_library
    Efield, Hfield = solve(solver, source, phantom)
    # Process results...
end
```

This provides significant performance benefits for parameter studies and optimization.

### Solver Extensions

The solver framework is designed for extensibility:

- **New source types**: Implement the `AbstractCurrentSource` interface
- **Additional solvers**: Follow the `AbstractElectromagneticSolver` pattern
- **Custom operators**: Extend the LinearSolve.jl operator interface
- **Performance backends**: Integrate with specialized linear algebra libraries

The modular design enables easy integration of new electromagnetic methods while maintaining the same user interface.