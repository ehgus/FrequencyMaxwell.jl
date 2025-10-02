"""
Medium type system for electromagnetic simulations.

This module provides abstract types and concrete implementations for representing
electromagnetic media with their permittivity distributions and background properties.
"""

"""
    AbstractMedium{T <: AbstractFloat}

Abstract base type for electromagnetic media.

All concrete medium types must implement:
- `permittivity(medium)`: Returns the permittivity distribution
- `permittivity_bg(medium)`: Returns the background permittivity
- `grid_size(medium)`: Returns the grid dimensions
"""
abstract type AbstractMedium{T <: AbstractFloat} end

"""
    Medium{T <: AbstractFloat, A <: AbstractArray{Complex{T}, 3}} <: AbstractMedium{T}

Concrete electromagnetic medium type with permittivity distribution.

This type encapsulates both the permittivity distribution and the background
permittivity value, providing a clean abstraction for electromagnetic simulations.

# Fields
- `permittivity::A`: 3D permittivity distribution (relative permittivity, εᵣ)
- `permittivity_bg::T`: Background permittivity (εᵣ_bg)

# Type Parameters
- `T`: Floating-point precision (Float32, Float64)
- `A`: Array type for permittivity (CPU/GPU flexibility)

# Example
```julia
# Create a medium with a spherical scatterer
grid_size = (128, 128, 64)
perm_distribution = ones(ComplexF64, grid_size)
# ... modify perm_distribution to add scatterer ...

medium = Medium(perm_distribution, 1.33^2)  # Water background

# Access properties
ε = permittivity(medium)
ε_bg = permittivity_bg(medium)
size = grid_size(medium)
```
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

Extract the permittivity distribution from a medium.
"""
permittivity(medium::Medium) = medium.permittivity

"""
    permittivity_bg(medium::AbstractMedium) -> T

Extract the background permittivity from a medium.
"""
permittivity_bg(medium::Medium) = medium.permittivity_bg

"""
    grid_size(medium::AbstractMedium) -> NTuple{3, Int}

Get the grid dimensions of a medium.
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
