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
abstract type AbstractElectromagneticSolver{T <: AbstractFloat} end

"""
    AbstractCurrentSource{T<:AbstractFloat}

Abstract supertype for all electromagnetic current sources.
The type parameter `T` specifies the floating-point precision used for calculations.
"""
abstract type AbstractCurrentSource{T <: AbstractFloat} end

