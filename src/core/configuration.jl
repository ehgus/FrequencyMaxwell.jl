"""
Configuration types for FrequencyMaxwell solvers.

This module defines immutable configuration structs that provide type-safe,
validated parameter storage for electromagnetic solvers.
"""

"""
    ConvergentBornConfig{T<:AbstractFloat} <: AbstractMaxwellConfig{T}

Immutable configuration for the Convergent Born Iterative Method solver.

This configuration struct replaces the mutable Dict-based approach from the legacy
implementation with a type-safe, validated parameter structure.

# Type Parameters
- `T`: Floating-point precision type (Float32 or Float64)

# Fields
- `wavelength::T`: Wavelength in the background medium (meters)
- `permittivity_bg::T`: Background relative permittivity
- `resolution::NTuple{3, T}`: Spatial resolution (dx, dy, dz) in meters  
- `grid_size::NTuple{3, Int}`: Number of grid points (Nx, Ny, Nz)
- `use_abbe_sine::Bool`: Whether to use Abbe sine condition for illumination
- `boundary_thickness::NTuple{3, T}`: PML boundary layer thickness in each direction
- `field_attenuation::NTuple{3, T}`: Field attenuation layer thickness in each direction
- `field_attenuation_sharpness::T`: Sharpness factor for field attenuation (0-1)
- `periodic_boundary::NTuple{3, Bool}`: Periodic boundary conditions (x, y, z)
- `iterations_max::Int`: Maximum number of Born iterations (-1 for auto)
- `tolerance::T`: Convergence tolerance for iterative solver
- `linear_solver::Symbol`: LinearSolve.jl algorithm (:iterative, :gmres, :bicgstab, :cg)
- `linear_solver_options::Dict{Symbol, Any}`: Additional solver options
- `preconditioner::Symbol`: Preconditioning strategy (:none, :diagonal, :ilu)

# Constructor
```julia
ConvergentBornConfig(;
    wavelength::Real,
    permittivity_bg::Real = 1.0,
    resolution::NTuple{3, <:Real},
    grid_size::NTuple{3, Int},
    use_abbe_sine::Bool = true,
    boundary_thickness::NTuple{3, <:Real} = (0.0, 0.0, 0.0),
    field_attenuation::NTuple{3, <:Real} = (0.0, 0.0, 0.0),
    field_attenuation_sharpness::Real = 1.0,
    periodic_boundary::NTuple{3, Bool} = (true, true, false),
    iterations_max::Int = -1,
    tolerance::Real = 1e-6
)
```

# Validation
The constructor automatically validates parameters:
- `wavelength > 0`
- `permittivity_bg > 0`
- `all(resolution .> 0)`
- `all(grid_size .> 0)`
- `0 ≤ tolerance ≤ 1`

# Example
```julia
config = ConvergentBornConfig(
    wavelength = 500e-9,     # 500 nm
    permittivity_bg = 1.33^2, # Water background
    resolution = (50e-9, 50e-9, 50e-9),  # 50 nm isotropic
    grid_size = (128, 128, 32)
)
```
"""
struct ConvergentBornConfig{T<:AbstractFloat} <: AbstractMaxwellConfig{T}
    wavelength::T
    permittivity_bg::T
    resolution::NTuple{3, T}
    grid_size::NTuple{3, Int}
    use_abbe_sine::Bool
    boundary_thickness::NTuple{3, T}
    field_attenuation::NTuple{3, T}
    field_attenuation_sharpness::T
    periodic_boundary::NTuple{3, Bool}
    iterations_max::Int
    tolerance::T
    linear_solver::Symbol
    linear_solver_options::Dict{Symbol, Any}
    preconditioner::Symbol
    
    function ConvergentBornConfig{T}(
        wavelength::T,
        permittivity_bg::T,
        resolution::NTuple{3, T},
        grid_size::NTuple{3, Int},
        use_abbe_sine::Bool,
        boundary_thickness::NTuple{3, T},
        field_attenuation::NTuple{3, T},
        field_attenuation_sharpness::T,
        periodic_boundary::NTuple{3, Bool},
        iterations_max::Int,
        tolerance::T,
        linear_solver::Symbol,
        linear_solver_options::Dict{Symbol, Any},
        preconditioner::Symbol
    ) where T<:AbstractFloat
        
        # Validate parameters
        wavelength > 0 || throw(ArgumentError("wavelength must be positive"))
        permittivity_bg > 0 || throw(ArgumentError("permittivity_bg must be positive"))
        all(resolution .> 0) || throw(ArgumentError("all resolution components must be positive"))
        all(grid_size .> 0) || throw(ArgumentError("all grid_size components must be positive"))
        0 ≤ tolerance ≤ 1 || throw(ArgumentError("tolerance must be in [0, 1]"))
        
        return new{T}(
            wavelength, permittivity_bg, resolution, grid_size,
            use_abbe_sine, boundary_thickness, field_attenuation, field_attenuation_sharpness,
            periodic_boundary, iterations_max, tolerance,
            linear_solver, linear_solver_options, preconditioner
        )
    end
end

# Main constructor with automatic type promotion
function ConvergentBornConfig(;
    wavelength::Real,
    permittivity_bg::Real = 1.0,
    resolution::NTuple{3, <:Real},
    grid_size::NTuple{3, Int},
    use_abbe_sine::Bool = true,
    boundary_thickness::NTuple{3, <:Real} = (0.0, 0.0, 0.0),
    field_attenuation::NTuple{3, <:Real} = (0.0, 0.0, 0.0),
    field_attenuation_sharpness::Real = 1.0,
    periodic_boundary::NTuple{3, Bool} = (true, true, false),
    iterations_max::Int = -1,
    tolerance::Real = 1e-6,
    linear_solver::Symbol = :iterative,
    linear_solver_options::Dict{Symbol, Any} = Dict{Symbol, Any}(),
    preconditioner::Symbol = :none
)
    # Promote to common floating-point type
    T = promote_type(typeof(wavelength), typeof(permittivity_bg), 
                     eltype(resolution), eltype(boundary_thickness), 
                     eltype(field_attenuation), typeof(field_attenuation_sharpness), typeof(tolerance))
    T <: AbstractFloat || (T = Float64)  # Fallback to Float64 if not floating-point
    
    return ConvergentBornConfig{T}(
        T(wavelength), T(permittivity_bg), 
        T.(resolution), grid_size, use_abbe_sine, 
        T.(boundary_thickness), T.(field_attenuation), T(field_attenuation_sharpness),
        periodic_boundary, iterations_max, T(tolerance),
        linear_solver, linear_solver_options, preconditioner
    )
end

"""
    grid_spacing(config::ConvergentBornConfig) -> NTuple{3, T}

Calculate the physical grid spacing in each direction.
"""
grid_spacing(config::ConvergentBornConfig) = config.resolution

"""
    domain_size(config::ConvergentBornConfig) -> NTuple{3, T}

Calculate the total physical domain size in each direction.
"""
function domain_size(config::ConvergentBornConfig{T}) where T
    return ntuple(i -> config.grid_size[i] * config.resolution[i], 3)
end

"""
    wavenumber_background(config::ConvergentBornConfig) -> T

Calculate the background medium wavenumber k₀ = 2π√(εᵦ)/λ₀.
"""
function wavenumber_background(config::ConvergentBornConfig{T}) where T
    return T(2π) * sqrt(config.permittivity_bg) / config.wavelength
end

"""
    show(io::IO, config::ConvergentBornConfig)

Custom display for configuration objects with formatted output.
"""
function Base.show(io::IO, config::ConvergentBornConfig{T}) where T
    print(io, "ConvergentBornConfig{$T}:")
    print(io, "\n  wavelength: $(config.wavelength)")
    print(io, "\n  permittivity_bg: $(config.permittivity_bg)")
    print(io, "\n  resolution: $(config.resolution)")
    print(io, "\n  grid_size: $(config.grid_size)")
    print(io, "\n  domain_size: $(domain_size(config))")
    print(io, "\n  k₀: $(wavenumber_background(config))")
end
