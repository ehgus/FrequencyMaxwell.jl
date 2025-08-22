"""
Convergent Born Iterative Method solver implementation.

This module provides the main electromagnetic solver using the Convergent Born
Iterative Method (CBS) with LinearSolve.jl integration and AD compatibility.
"""

"""
    ConvergentBornSolver{T, AT} <: AbstractElectromagneticSolver{T}

Mutable solver state for the Convergent Born Iterative Method.

This struct separates immutable configuration from mutable solver state,
following modern Julia best practices for performance and AD compatibility.

# Type Parameters
- `T<:AbstractFloat`: Floating-point precision type
- `AT<:AbstractArray`: Array type (for CPU/GPU flexibility)

# Fields
- `config::ConvergentBornConfig{T}`: Immutable solver configuration
- `Green_function::Union{Nothing, AbstractArray}`: Cached Green's function
- `potential::Union{Nothing, AT}`: Current iteration potential (V = δε)
- `internal_fields::Union{Nothing, AT}`: Cached internal electromagnetic fields
- `iteration_count::Int`: Current iteration number
- `residual_history::Vector{T}`: Convergence history tracking
- `solver_state::Any`: LinearSolve.jl solver state for performance

# Constructor
```julia
ConvergentBornSolver(config::ConvergentBornConfig{T}; array_type=Array) where T
```

# Example
```julia
config = ConvergentBornConfig(
    wavelength = 500e-9,
    NA = 1.4,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32)
)

solver = ConvergentBornSolver(config)
```
"""
mutable struct ConvergentBornSolver{T<:AbstractFloat, AT<:AbstractArray} <: AbstractElectromagneticSolver{T}
    config::ConvergentBornConfig{T}
    Green_function::Union{Nothing, AbstractArray{Complex{T}}}
    potential::Union{Nothing, AT}
    internal_fields::Union{Nothing, AT}
    iteration_count::Int
    residual_history::Vector{T}
    solver_state::Any  # LinearSolve.jl solver state
    
    function ConvergentBornSolver{T, AT}(
        config::ConvergentBornConfig{T}
    ) where {T<:AbstractFloat, AT<:AbstractArray}
        new{T, AT}(
            config,
            nothing,  # Green_function - computed lazily
            nothing,  # potential - set during solve
            nothing,  # internal_fields - computed during solve
            0,        # iteration_count
            T[],      # residual_history
            nothing   # solver_state - initialized during solve
        )
    end
end

# Convenience constructor with automatic array type inference
function ConvergentBornSolver(
    config::ConvergentBornConfig{T};
    array_type::Type{<:AbstractArray} = Array
) where T<:AbstractFloat
    AT = array_type{Complex{T}, 3}  # 3D complex arrays for electromagnetic fields
    return ConvergentBornSolver{T, AT}(config)
end

"""
    solve(solver::ConvergentBornSolver, source::AbstractCurrentSource, permittivity::AbstractArray)

Solve the electromagnetic scattering problem using the Convergent Born method.

This function implements the core CBS iteration with LinearSolve.jl integration
for modern, flexible linear algebra backends and AD compatibility.

# Arguments
- `solver::ConvergentBornSolver`: Pre-configured solver instance
- `source::AbstractCurrentSource`: Electromagnetic current source
- `permittivity::AbstractArray{T, 3}`: 3D permittivity distribution

# Returns
- `E_field::AbstractArray{Complex{T}, 4}`: Electric field (last dim = 3 for components)
- `H_field::AbstractArray{Complex{T}, 4}`: Magnetic field (last dim = 3 for components)

# Algorithm
The method solves the electromagnetic scattering equation:
```
ψ = ψ_incident + G * V * ψ
```
where:
- ψ: total field
- ψ_incident: incident field from source
- G: dyadic Green's function
- V: potential (permittivity contrast)

# LinearSolve.jl Integration
The iteration is formulated as a linear system:
```
(I - G*V) * ψ = ψ_incident
```
This enables use of modern iterative solvers and preconditioning.

# Example
```julia
# Create permittivity distribution
perm = ones(Complex{Float64}, solver.config.grid_size...)
perm[50:80, 50:80, 10:20] .= 1.5^2  # Add scatterer

# Solve
E_field, H_field = solve(solver, source, perm)
```
"""
function solve(
    solver::ConvergentBornSolver{T, AT},
    source::AbstractCurrentSource{T},
    permittivity::AbstractArray{<:Number, 3}
) where {T<:AbstractFloat, AT<:AbstractArray}
    
    # Validate input dimensions
    size(permittivity) == solver.config.grid_size || 
        throw(ArgumentError("permittivity size $(size(permittivity)) must match grid_size $(solver.config.grid_size)"))
    
    # Initialize solver state
    _initialize_solver!(solver, permittivity)
    
    # Generate incident fields from source
    E_incident, H_incident = _generate_incident_fields(solver, source)
    
    # Solve scattering equation using LinearSolve.jl
    E_total, H_total = _solve_scattering_equation(solver, E_incident, H_incident)
    
    # Update solver statistics
    solver.iteration_count += 1
    
    return E_total, H_total
end

"""
    _initialize_solver!(solver::ConvergentBornSolver, permittivity::AbstractArray)

Initialize solver internal state for a new problem.
"""
function _initialize_solver!(
    solver::ConvergentBornSolver{T, AT},
    permittivity::AbstractArray{<:Number, 3}
) where {T<:AbstractFloat, AT<:AbstractArray}
    
    # Compute potential (permittivity contrast)
    solver.potential = AT(Complex{T}.(permittivity .- solver.config.permittivity_bg))
    
    # Initialize Green's function if needed
    if solver.Green_function === nothing
        solver.Green_function = _compute_green_function(solver)
    end
    
    # Reset iteration tracking
    solver.iteration_count = 0
    empty!(solver.residual_history)
    
    return nothing
end

"""
    _compute_green_function(solver::ConvergentBornSolver) -> AbstractArray

Compute the dyadic Green's function for the given configuration.
This is a placeholder for the actual Green's function computation.
"""
function _compute_green_function(solver::ConvergentBornSolver{T}) where T
    # Placeholder implementation - actual implementation would compute
    # the electromagnetic Green's function based on solver configuration
    
    grid_size = solver.config.grid_size
    # Return dummy Green's function array for now
    return zeros(Complex{T}, grid_size..., 3, 3)  # 3x3 dyadic for each grid point
end

"""
    _generate_incident_fields(solver::ConvergentBornSolver, source::AbstractCurrentSource)

Generate incident electromagnetic fields from the current source.
"""
function _generate_incident_fields(
    solver::ConvergentBornSolver{T, AT},
    source::AbstractCurrentSource{T}
) where {T<:AbstractFloat, AT<:AbstractArray}
    
    # Placeholder implementation - actual implementation would compute
    # incident fields based on source type and solver configuration
    
    grid_size = solver.config.grid_size
    E_incident = zeros(Complex{T}, grid_size..., 3)
    H_incident = zeros(Complex{T}, grid_size..., 3)
    
    # Simple plane wave implementation as placeholder
    if isa(source, PlaneWaveSource)
        _fill_plane_wave_fields!(E_incident, H_incident, solver, source)
    end
    
    return E_incident, H_incident
end

"""
    _fill_plane_wave_fields!(E_incident, H_incident, solver, source)

Fill incident field arrays with plane wave solution.
"""
function _fill_plane_wave_fields!(
    E_incident::AbstractArray{Complex{T}, 4},
    H_incident::AbstractArray{Complex{T}, 4},
    solver::ConvergentBornSolver{T},
    source::PlaneWaveSource{T}
) where T<:AbstractFloat
    
    # Placeholder - actual implementation would compute plane wave fields
    # based on source parameters and grid configuration
    
    # For now, set simple test field
    E_incident[:, :, :, 1] .= Complex{T}(1.0)  # Ex = 1
    H_incident[:, :, :, 2] .= Complex{T}(1.0/377.0)  # Hy = Ex/Z0
    
    return nothing
end

"""
    _solve_scattering_equation(solver, E_incident, H_incident)

Solve the electromagnetic scattering equation using LinearSolve.jl.
"""
function _solve_scattering_equation(
    solver::ConvergentBornSolver{T, AT},
    E_incident::AbstractArray{Complex{T}, 4},
    H_incident::AbstractArray{Complex{T}, 4}
) where {T<:AbstractFloat, AT<:AbstractArray}
    
    # For now, return incident fields as placeholder
    # Actual implementation would:
    # 1. Set up LinearProblem with (I - G*V) system matrix
    # 2. Use LinearSolve.jl to solve the system
    # 3. Return total fields
    
    return E_incident, H_incident
end

"""
    reset!(solver::ConvergentBornSolver)

Reset solver state for a new problem while preserving configuration.
"""
function reset!(solver::ConvergentBornSolver)
    solver.potential = nothing
    solver.internal_fields = nothing
    solver.iteration_count = 0
    empty!(solver.residual_history)
    solver.solver_state = nothing
    return nothing
end

"""
    show(io::IO, solver::ConvergentBornSolver)

Custom display for solver objects with status information.
"""
function Base.show(io::IO, solver::ConvergentBornSolver{T, AT}) where {T, AT}
    print(io, "ConvergentBornSolver{$T, $(AT.name)}:")
    print(io, "\n  iterations: $(solver.iteration_count)")
    if !isempty(solver.residual_history)
        print(io, "\n  last residual: $(last(solver.residual_history))")
    end
    print(io, "\n  configuration:")
    show(io, solver.config)
end
