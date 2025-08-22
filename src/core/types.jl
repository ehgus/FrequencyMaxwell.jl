"""
Core type definitions for FrequencyMaxwell.

This module defines the fundamental abstract types and type aliases used throughout
the package, following modern Julia type system best practices.
"""

"""
    AbstractElectromagneticSolver{T<:AbstractFloat}

Abstract supertype for all electromagnetic solvers in FrequencyMaxwell.
The type parameter `T` specifies the floating-point precision used for calculations.
"""
abstract type AbstractElectromagneticSolver{T<:AbstractFloat} end

"""
    AbstractCurrentSource{T<:AbstractFloat}

Abstract supertype for all electromagnetic current sources.
The type parameter `T` specifies the floating-point precision used for calculations.
"""
abstract type AbstractCurrentSource{T<:AbstractFloat} end

"""
    AbstractMaxwellConfig{T<:AbstractFloat}

Abstract supertype for solver configuration objects.
The type parameter `T` specifies the floating-point precision used for calculations.
"""
abstract type AbstractMaxwellConfig{T<:AbstractFloat} end

"""
Type aliases for common array types to ensure consistency across the package.
"""
const RealVec3{T} = SVector{3, T} where T<:Real
const ComplexVec3{T} = SVector{3, Complex{T}} where T<:Real
const RealArray3D{T} = AbstractArray{T, 3} where T<:Real
const ComplexArray3D{T} = AbstractArray{Complex{T}, 3} where T<:Real
const ComplexArray4D{T} = AbstractArray{Complex{T}, 4} where T<:Real
