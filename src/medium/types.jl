"""
Medium type system for electromagnetic simulations.

This module provides abstract types and concrete implementations for representing
electromagnetic media with their permittivity distributions and background properties.

# Design Rationale

The Medium type abstracts material properties in electromagnetic simulations, providing:

1. **Type Safety**: Parametric types ensure compile-time verification of precision and array types
2. **GPU Compatibility**: Array type parameter `A` supports both CPU and GPU arrays
3. **Immutability**: Struct fields are immutable for thread safety, though array contents remain mutable
4. **Validation**: Constructor validates physical constraints (positive permittivity)
5. **Performance**: Type parameters enable compiler optimizations and zero-cost abstractions

# Mathematical Foundation

The Medium type represents the relative permittivity distribution εᵣ(r) and background
permittivity εᵣ_bg, which define the electromagnetic scattering potential:

```math
V(r) = k₀² [εᵣ(r)/εᵣ_bg - 1]
```

where k₀ = 2π/λ₀ is the free-space wavenumber.

See also: [`AbstractMedium`](@ref), [`Medium`](@ref), [`permittivity`](@ref),
[`permittivity_bg`](@ref), [`grid_size`](@ref)
"""

"""
    AbstractMedium{T <: AbstractFloat}

Abstract base type for electromagnetic media in frequency-domain simulations.

# Type Parameter

- `T`: Floating-point precision (Float32 or Float64)

# Interface Requirements

All concrete subtypes must implement:
- `permittivity(medium)`: Returns the 3D permittivity distribution εᵣ(r)
- `permittivity_bg(medium)`: Returns the background permittivity εᵣ_bg
- `grid_size(medium)`: Returns the grid dimensions (Nx, Ny, Nz)

# Design Rationale

The abstract type enables:
1. Multiple medium representations (homogeneous, inhomogeneous, anisotropic)
2. Consistent interface for electromagnetic solvers
3. Type-stable dispatch based on precision level
4. Future extensibility for advanced material models

See also: [`Medium`](@ref)
"""
abstract type AbstractMedium{T <: AbstractFloat} end

"""
    Medium{T <: AbstractFloat, A <: AbstractArray{Complex{T}, 3}} <: AbstractMedium{T}

Concrete electromagnetic medium type with 3D permittivity distribution.

Represents an electromagnetic medium characterized by its relative permittivity
distribution εᵣ(r) on a regular 3D grid, along with a background permittivity εᵣ_bg.

# Fields

- `permittivity::A`: 3D complex permittivity distribution εᵣ(r) = εᵣ' + iεᵣ''
  - Real part (εᵣ') represents dielectric properties
  - Imaginary part (εᵣ'') represents absorption/loss
- `permittivity_bg::T`: Background permittivity εᵣ_bg (must be positive)

# Type Parameters

- `T`: Floating-point precision (Float32 or Float64)
- `A`: Array type (supports CPU arrays, CuArray, ROCArray, MtlArray, etc.)

# Constructor

    Medium(permittivity::AbstractArray{Complex{T}, 3}, permittivity_bg::Real)

Validates that `permittivity_bg > 0` and constructs a Medium instance.

# Physical Interpretation

The relative permittivity relates to refractive index via:
```math
n(r) = \\sqrt{εᵣ(r)}
```

For lossy media:
```math
εᵣ(r) = εᵣ'(r) + iεᵣ''(r)
```
where εᵣ'' > 0 indicates absorption.

# Performance Considerations

1. **Array Type Flexibility**: Type parameter `A` enables GPU acceleration
   - CPU: `Array{ComplexF64, 3}`
   - CUDA: `CuArray{ComplexF64, 3}`
   - ROCm: `ROCArray{ComplexF64, 3}`
   - Metal: `MtlArray{ComplexF64, 3}`

2. **Type Stability**: Parametric types enable compiler optimizations

3. **Memory Layout**: Permittivity arrays should use column-major (Julia native) layout

# Examples

```julia
# Basic usage: homogeneous medium
using FrequencyMaxwell

gs = (64, 64, 32)
perm = ones(ComplexF64, gs) .* 2.25  # Glass (n ≈ 1.5)
medium = Medium(perm, 1.0)  # Air background

# Inhomogeneous medium: spherical scatterer
perm = ones(ComplexF64, 128, 128, 64)
# Add sphere at center with higher permittivity
center = (64, 64, 32)
radius = 10
for i in 1:128, j in 1:128, k in 1:64
    if sqrt((i-center[1])^2 + (j-center[2])^2 + (k-center[3])^2) < radius
        perm[i, j, k] = 4.0 + 0.0im  # n = 2.0
    end
end
medium = Medium(perm, 1.33^2)  # Water background

# Lossy medium with absorption
perm_lossy = (2.25 + 0.1im) .* ones(ComplexF64, 64, 64, 32)
medium_lossy = Medium(perm_lossy, 1.0)

# GPU-accelerated medium (requires CUDA.jl)
# using CUDA
# perm_gpu = CUDA.CuArray(perm)
# medium_gpu = Medium(perm_gpu, 1.0)

# Mixed precision for memory efficiency
perm_f32 = ones(ComplexF32, 128, 128, 64)
medium_f32 = Medium(perm_f32, 1.0f0)  # Uses Float32

# Access medium properties
ε = permittivity(medium)        # Get permittivity distribution
ε_bg = permittivity_bg(medium)  # Get background permittivity
size = grid_size(medium)        # Get grid dimensions (Nx, Ny, Nz)
```

# See Also

- [`AbstractMedium`](@ref): Abstract base type
- [`permittivity`](@ref): Extract permittivity distribution
- [`permittivity_bg`](@ref): Extract background permittivity
- [`grid_size`](@ref): Extract grid dimensions
- [`ConvergentBornSolver`](@ref): Solver that uses Medium types
"""
struct Medium{T <: AbstractFloat, A <: AbstractArray{Complex{T}, 3}} <: AbstractMedium{T}
    permittivity::A
    permittivity_bg::T

    function Medium(
            permittivity::A,
            permittivity_bg::Real
    ) where {T <: AbstractFloat, A <: AbstractArray{Complex{T}, 3}}
        # Validate inputs
        permittivity_bg > 0 ||
            throw(ArgumentError("permittivity_bg must be positive, got: $permittivity_bg"))

        new{T, A}(permittivity, T(permittivity_bg))
    end
end

# Convenience constructor for automatic type inference
function Medium(
        permittivity::AbstractArray{<:Complex, 3},
        permittivity_bg::Real
)
    T = real(eltype(permittivity))
    return Medium{T, typeof(permittivity)}(permittivity, permittivity_bg)
end

"""
    permittivity(medium::AbstractMedium) -> AbstractArray{Complex{T}, 3}

Extract the 3D permittivity distribution εᵣ(r) from a medium.

Returns the complex relative permittivity array where:
- Real part represents dielectric properties (εᵣ' ≥ 1 typically)
- Imaginary part represents absorption/loss (εᵣ'' ≥ 0)

# Example

```julia
medium = Medium(ones(ComplexF64, 64, 64, 32), 1.0)
ε = permittivity(medium)
ε_real = real(ε)  # Dielectric component
ε_imag = imag(ε)  # Loss component
```

See also: [`permittivity_bg`](@ref), [`Medium`](@ref)
"""
permittivity(medium::Medium) = medium.permittivity

"""
    permittivity_bg(medium::AbstractMedium) -> T

Extract the background permittivity εᵣ_bg from a medium.

Returns the real-valued background relative permittivity, which represents the
medium surrounding the scattering region. This value is always positive and
defines the reference refractive index n_bg = √εᵣ_bg.

# Example

```julia
medium = Medium(ones(ComplexF64, 64, 64, 32), 1.33^2)  # Water background
ε_bg = permittivity_bg(medium)  # Returns 1.7689
n_bg = sqrt(ε_bg)                # Returns 1.33 (water refractive index)
```

See also: [`permittivity`](@ref), [`Medium`](@ref)
"""
permittivity_bg(medium::Medium) = medium.permittivity_bg

"""
    grid_size(medium::AbstractMedium) -> NTuple{3, Int}

Get the grid dimensions (Nx, Ny, Nz) of a medium's permittivity distribution.

Returns a tuple of three integers representing the number of grid points in
the x, y, and z directions.

# Example

```julia
medium = Medium(ones(ComplexF64, 128, 64, 32), 1.0)
(Nx, Ny, Nz) = grid_size(medium)  # Returns (128, 64, 32)
```

See also: [`Medium`](@ref), [`permittivity`](@ref)
"""
grid_size(medium::Medium) = size(medium.permittivity)

"""
    Base.show(io::IO, medium::Medium{T}) where T

Display method for Medium type.
"""
function Base.show(io::IO, medium::Medium{T}) where {T}
    gs = grid_size(medium)
    perm_bg = permittivity_bg(medium)
    n_bg = sqrt(perm_bg)

    # Calculate permittivity range
    perm_values = medium.permittivity
    perm_min = minimum(abs.(perm_values))
    perm_max = maximum(abs.(perm_values))
    n_min = sqrt(perm_min)
    n_max = sqrt(perm_max)

    println(io, "Medium{$T}")
    println(io, "  Grid size: $(gs[1]) × $(gs[2]) × $(gs[3])")
    println(io, "  Background: εᵣ = $(round(perm_bg, digits=4)), n = $(round(n_bg, digits=4))")
    println(io, "  Permittivity range: $(round(perm_min, digits=4)) to $(round(perm_max, digits=4))")
    print(io, "  Refractive index range: $(round(n_min, digits=4)) to $(round(n_max, digits=4))")
end

# Export types and functions
export AbstractMedium, Medium
export permittivity, permittivity_bg, grid_size
