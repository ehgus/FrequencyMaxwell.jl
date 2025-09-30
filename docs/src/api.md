# API Reference

```@meta
CurrentModule = FrequencyMaxwell
```

This page documents the complete API for FrequencyMaxwell.jl, organized by functionality.

## Solver Types

### ConvergentBornSolver

The main electromagnetic solver with integrated configuration.

```@docs
ConvergentBornSolver
```

### Configuration Utilities

Utility functions for working with solver parameters.

```@docs
domain_size
grid_spacing
wavenumber_background
```

### Core Solving Function

```@docs
solve
```

## Source Types

### Abstract Base Types

```@docs
AbstractCurrentSource
```

### Concrete Source Implementations

```@docs
PlaneWaveSource
```

### Source Utilities

```@docs
source_wavelength
source_power
validate_source
generate_incident_fields
```

## Field Utilities

### ElectromagneticField Type

```@docs
ElectromagneticField
```

### Field Analysis Functions

```@docs
field_energy
field_intensity
poynting_vector
extract_plane
```

## Geometry and Phantoms

### Phantom Generation Functions

```@docs
phantom_bead
phantom_plate
phantom_cylinder
```

## Internal Core Types

These types are primarily for internal use but may be useful for advanced users extending the package.

### Abstract Base Types

The following abstract types form the foundation of the type hierarchy:

- `AbstractElectromagneticSolver{T<:AbstractFloat}`: Base type for all electromagnetic solvers
- `AbstractMaxwellConfig{T<:AbstractFloat}`: Base type for all configuration objects

### Type Aliases

The package defines several type aliases for consistency:

- `RealVec3{T}`: 3D real vector using StaticArrays
- `ComplexVec3{T}`: 3D complex vector using StaticArrays
- `RealArray3D{T}`: 3D real array (hardware-agnostic)
- `ComplexArray3D{T}`: 3D complex array (hardware-agnostic)
- `ComplexArray4D{T}`: 4D complex array (hardware-agnostic)

### Green's Function Implementation

The core electromagnetic kernels are implemented in:

- `DyadicGreen{T}`: Green's function operator for electromagnetic propagation
- `conv!`: In-place convolution operation using FFT

These components handle the mathematical foundation of the Convergent Born Iterative Method (CBS) and provide hardware-agnostic computation through KernelAbstractions.jl.

## LinearSolve.jl Integration

FrequencyMaxwell.jl integrates deeply with LinearSolve.jl for flexible linear algebra backends:

### Custom Linear Operators

- `CBSLinearOperator`: Implements the (I - G*V) operator for CBS reformulation

### Supported Algorithms

The package supports any LinearSolve.jl algorithm, with tested compatibility for:
- `KrylovJL_GMRES()`: Recommended default Krylov method
- `GMRES()`: Generalized Minimal RESidual method
- `BiCGStabL()`: Bi-Conjugate Gradient Stabilized
- `IDRsAlgorithm()`: Induced Dimension Reduction
- And many others from the LinearSolve.jl ecosystem

## Index

```@index
```