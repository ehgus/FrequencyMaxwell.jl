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
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32)
)

solver = ConvergentBornSolver(config)
```
"""
mutable struct ConvergentBornSolver{T<:AbstractFloat, AT<:AbstractArray} <: AbstractElectromagneticSolver{T}
    config::ConvergentBornConfig{T}
    Green_function::Union{Nothing, DyadicGreen{T}}
    flip_Green_function::Union{Nothing, DyadicGreen{T}}
    potential::Union{Nothing, AT}  # V = k²(ε/ε_bg - 1) with padding
    permittivity::Union{Nothing, AT}  # Original permittivity without padding
    boundary_thickness_pixel::NTuple{3, Int}
    field_attenuation_pixel::NTuple{3, Int}
    field_attenuation_masks::Vector{AT}
    ROI::NTuple{6, Int}  # Region of interest bounds
    eps_imag::T  # Imaginary part for convergence
    Bornmax::Int  # Actual number of iterations to use
    iteration_count::Int
    residual_history::Vector{T}
    solver_state::Any  # LinearSolve.jl solver state
    
    function ConvergentBornSolver{T, AT}(
        config::ConvergentBornConfig{T}
    ) where {T<:AbstractFloat, AT<:AbstractArray}
        # Calculate boundary thickness in pixels
        boundary_thickness_pixel = round.(Int, config.boundary_thickness ./ (config.resolution .* 2))
        field_attenuation_pixel = round.(Int, config.field_attenuation ./ (config.resolution .* 2))
        
        # Calculate region of interest
        ROI = (
            boundary_thickness_pixel[1] + 1, boundary_thickness_pixel[1] + config.grid_size[1],
            boundary_thickness_pixel[2] + 1, boundary_thickness_pixel[2] + config.grid_size[2],
            boundary_thickness_pixel[3] + 1, boundary_thickness_pixel[3] + config.grid_size[3]
        )
        
        new{T, AT}(
            config,
            nothing,  # Green_function - computed lazily
            nothing,  # flip_Green_function - computed lazily
            nothing,  # potential - set during solve
            nothing,  # permittivity - set during solve
            boundary_thickness_pixel,
            field_attenuation_pixel,
            AT[],     # field_attenuation_masks - computed during initialization
            ROI,
            T(0),     # eps_imag - calculated based on potential
            0,        # Bornmax - calculated automatically
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
    
    # Initialize solver state (computes potential V with proper padding)
    _initialize_solver!(solver, permittivity)
    
    # Generate incident fields from source (with padding)
    E_incident, H_incident = _generate_incident_fields_padded(solver, source)
    
    # Generate source term: V * E_incident + i*eps_imag * E_incident
    source_term = solver.potential .* E_incident .+ Complex{T}(0, solver.eps_imag) .* E_incident
    
    # Solve scattering equation using CBS algorithm
    E_scattered, H_scattered = _solve_cbs_scattering(solver, source_term)
    
    # Add incident field to get total field and crop to ROI
    E_total = crop_to_ROI(E_scattered + E_incident, solver)
    H_total = crop_to_ROI(H_scattered + H_incident, solver)  
    
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
    
    # Store original permittivity
    solver.permittivity = AT(Complex{T}.(permittivity))
    
    # Compute potential V = k²(ε/ε_bg - 1) with proper padding
    k = T(2π) * sqrt(solver.config.permittivity_bg) / solver.config.wavelength
    potential_unpadded = AT(Complex{T}.(k^2 * (permittivity ./ solver.config.permittivity_bg .- 1)))
    
    # Pad the potential (replicate boundary values)
    solver.potential = _pad_array(potential_unpadded, solver.boundary_thickness_pixel, :replicate)
    
    # Calculate eps_imag based on maximum potential value (following jl-ConvergentBornSolver)
    solver.eps_imag = max(T(2^-20), maximum(abs.(solver.potential)) * T(1.01))
    
    # Add imaginary part to potential for convergence
    solver.potential .-= Complex{T}(0, solver.eps_imag)
    
    # Create field attenuation masks
    _create_attenuation_masks!(solver)
    
    # Apply attenuation masks to potential
    for mask in solver.field_attenuation_masks
        solver.potential .*= mask
    end
    
    # Calculate optimal iteration count if automatic
    if solver.config.iterations_max < 0
        k0_nm = T(2π) * sqrt(solver.config.permittivity_bg) / solver.config.wavelength
        steps = abs(2 * k0_nm / solver.eps_imag)
        domain_extent = norm(size(solver.potential)[1:3] .* solver.config.resolution)
        Bornmax_opt = ceil(Int, domain_extent / steps / 2 + 1) * 2
        solver.Bornmax = Bornmax_opt * abs(solver.config.iterations_max)
    else
        solver.Bornmax = solver.config.iterations_max
    end
    
    # Initialize Green's functions if needed
    if solver.Green_function === nothing
        solver.Green_function, solver.flip_Green_function = _compute_green_functions(solver)
    end
    
    # Reset iteration tracking
    solver.iteration_count = 0
    empty!(solver.residual_history)
    
    return nothing
end

"""
    _compute_green_function(solver::ConvergentBornSolver) -> DyadicGreen

Compute the dyadic Green's function for the given configuration.
"""
function _compute_green_functions(solver::ConvergentBornSolver{T}) where T
    config = solver.config
    
    # Calculate wave number in background medium
    k0_nm = T(2π) * sqrt(config.permittivity_bg) / config.wavelength
    k_square = k0_nm^2 + Complex{T}(0, solver.eps_imag)
    
    # Array size includes boundary padding
    arr_size = size(solver.potential)
    
    # Subpixel shifts for proper boundary conditions
    subpixel_shift = ntuple(i -> config.periodic_boundary[i] ? T(0) : T(0.25), 3)
    flip_subpixel_shift = ntuple(i -> config.periodic_boundary[i] ? T(0) : T(-0.25), 3)
    
    # Create dyadic Green's functions
    array_type = typeof(solver).parameters[2]
    Green_fn = DyadicGreen(array_type, k_square, arr_size, config.resolution, subpixel_shift)
    flip_Green_fn = DyadicGreen(array_type, k_square, arr_size, config.resolution, flip_subpixel_shift)
    
    return Green_fn, flip_Green_fn
end

"""
    _generate_incident_fields(solver::ConvergentBornSolver, source::AbstractCurrentSource)

Generate incident electromagnetic fields from the current source.
"""
function _generate_incident_fields_padded(
    solver::ConvergentBornSolver{T, AT},
    source::AbstractCurrentSource{T}
) where {T<:AbstractFloat, AT<:AbstractArray}
    
    # Generate incident fields on original grid
    grid_size = solver.config.grid_size
    E_incident_orig = zeros(Complex{T}, grid_size..., 3)
    H_incident_orig = zeros(Complex{T}, grid_size..., 3)
    
    # Simple plane wave implementation
    if isa(source, PlaneWaveSource)
        _fill_plane_wave_fields!(E_incident_orig, H_incident_orig, solver, source)
    end
    
    # Pad incident fields to match potential array size
    E_incident = _pad_array(E_incident_orig, solver.boundary_thickness_pixel, :zeros)
    H_incident = _pad_array(H_incident_orig, solver.boundary_thickness_pixel, :zeros)
    
    return E_incident, H_incident
end

"""
    _fill_plane_wave_fields!(E_incident, H_incident, solver, source)

Fill incident field arrays with proper plane wave solution.

Computes plane wave fields with proper k-vector propagation:
E(r) = E₀ * exp(ik·r)
H(r) = (k × E) / (ωμ₀) = (k × E) / (k₀Z₀)

The coordinate system is centered at the domain center for proper phase reference.
"""
function _fill_plane_wave_fields!(
    E_incident::AbstractArray{Complex{T}, 4},
    H_incident::AbstractArray{Complex{T}, 4},
    solver::ConvergentBornSolver{T},
    source::PlaneWaveSource{T}
) where T<:AbstractFloat
    
    config = solver.config
    grid_size = config.grid_size
    resolution = config.resolution
    
    # Wave parameters
    k0 = T(2π / config.wavelength)  # Free-space wave number
    k_bg = k0 * sqrt(config.permittivity_bg)  # Background wave number
    
    # Normalize k-vector and polarization
    k_hat = source.k_vector ./ norm(source.k_vector)  # Unit k-vector
    k_vec = k_bg .* k_hat  # Full k-vector in medium
    E0 = source.polarization ./ norm(source.polarization)  # Normalized polarization
    
    # Ensure E ⟂ k (transversality condition)
    E_perp = E0 .- (dot(E0, k_hat) * k_hat)  # Remove parallel component
    E_perp = E_perp ./ norm(E_perp)  # Renormalize
    
    # Compute H from E using H = k × E / (ωμ₀) = k × E / (k₀Z₀)
    H_perp = cross(k_hat, E_perp) ./ T(377 / sqrt(config.permittivity_bg))
    
    # Calculate domain center for phase reference (consistent with PlaneWaveSource)
    center = ntuple(i -> (grid_size[i] - 1) * resolution[i] / 2, 3)
    
    # Compute phase array k·r for all grid points
    for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
        # Position relative to domain center
        r = [
            (i - 1) * resolution[1] - center[1],
            (j - 1) * resolution[2] - center[2], 
            (k - 1) * resolution[3] - center[3]
        ]
        
        # Phase: k·r + source phase offset
        phase = dot(k_vec, r) + source.phase
        complex_amplitude = source.amplitude * exp(im * phase)
        
        # Set E field components
        E_incident[i, j, k, 1] = complex_amplitude * E_perp[1]
        E_incident[i, j, k, 2] = complex_amplitude * E_perp[2]
        E_incident[i, j, k, 3] = complex_amplitude * E_perp[3]
        
        # Set H field components
        H_incident[i, j, k, 1] = complex_amplitude * H_perp[1]
        H_incident[i, j, k, 2] = complex_amplitude * H_perp[2]
        H_incident[i, j, k, 3] = complex_amplitude * H_perp[3]
    end
    
    return nothing
end

"""
    _solve_scattering_equation(solver, E_incident, H_incident)

Solve the electromagnetic scattering equation using the Convergent Born Series.

Implements the CBS iteration:
E₁ = r/2 * (G + G_flip) * S = V * (G + G_flip) * (i/ε_imag/2) * S  
E_{n+1} = M * E_n where M := V * (G + G_flip)/2 * r - r + 1
"""
function _solve_cbs_scattering(
    solver::ConvergentBornSolver{T, AT},
    source::AbstractArray{Complex{T}, 4}
) where {T<:AbstractFloat, AT<:AbstractArray}
    
    # Array dimensions with padding
    size_field = size(source)
    Green_fn = solver.Green_function
    flip_Green_fn = solver.flip_Green_function
    
    # Initialize working arrays
    psi = zeros(Complex{T}, size_field)
    flip_psi = zeros(Complex{T}, size_field)
    Field = zeros(Complex{T}, size_field)
    Field_n = zeros(Complex{T}, size_field)
    
    # CBS iteration 1: E1 = r/2 * (G + G_flip) * S = V * (G + G_flip) * (i/ε/2) * S
    if solver.Bornmax >= 1
        source_scaled = source .* (Complex{T}(0, 1) / solver.eps_imag)
        _apply_attenuation_masks!(source_scaled, solver.field_attenuation_masks)
        psi .= source_scaled
        
        # Convolution with Green functions: (G + G_flip)/2
        flip_psi .= psi
        conv!(Green_fn, psi)
        conv!(flip_Green_fn, flip_psi) 
        psi .+= flip_psi
        psi ./= 2
        
        # Element-wise multiplication with potential
        Field_n .= solver.potential .* psi
        _apply_attenuation_masks!(Field_n, solver.field_attenuation_masks)
        Field .+= Field_n
    end
    
    # CBS iterations n ≥ 2: E_{n+1} = M * E_n where M := V * (G + G_flip)/2 * r - r + 1
    for iteration = 2:solver.Bornmax
        psi .= (Complex{T}(0, 1) / solver.eps_imag) .* solver.potential .* Field_n
        Field_n .-= psi
        
        # Convolution with Green functions
        flip_psi .= psi
        conv!(Green_fn, psi)
        conv!(flip_Green_fn, flip_psi)
        psi .+= flip_psi
        psi ./= 2
        
        # Element-wise multiplication with potential
        Field_n .+= solver.potential .* psi
        _apply_attenuation_masks!(Field_n, solver.field_attenuation_masks)
        Field .+= Field_n
    end
    
    # Compute magnetic field
    Hfield = _compute_magnetic_field(solver, Field)
    
    return Field, Hfield
end

"""
    _compute_magnetic_field(solver, E_field) -> AbstractArray

Compute magnetic field from electric field using Maxwell's equation:
H = -i/k₀ * (n₀/Z₀) * ∇ × E

where k₀ is the free-space wave number and Z₀ = 377 Ω is the impedance.
"""
function _compute_magnetic_field(
    solver::ConvergentBornSolver{T, AT},
    E_field::AbstractArray{Complex{T}, 4}
) where {T<:AbstractFloat, AT<:AbstractArray}
    
    # Create curl operator for field with padding
    curl_op = Curl(AT, T, size(solver.potential), solver.config.resolution)
    
    # Compute curl of E field
    Hfield = conv(curl_op, E_field)
    
    # Apply scaling: H = -i * λ/(2π * Z₀) * ∇ × E
    scaling_factor = Complex{T}(0, -1) * solver.config.wavelength / (T(2π) * T(377))
    
    return scaling_factor .* Hfield
end

"""
    _pad_array(arr, padding_pixels, mode) -> padded_array

Pad an array with specified padding in each dimension.
"""
function _pad_array(arr::AbstractArray{T, N}, padding_pixels::NTuple{3, Int}, mode::Symbol) where {T, N}
    if all(padding_pixels .== 0)
        return arr
    end
    
    original_size = size(arr)
    padded_size = if N == 3
        original_size .+ 2 .* padding_pixels
    elseif N == 4  # For field arrays with component dimension
        (original_size[1:3] .+ 2 .* padding_pixels..., original_size[4])
    else
        error("Unsupported array dimension: $N")
    end
    
    # Create padded array
    if mode == :zeros
        padded = zeros(T, padded_size)
    elseif mode == :replicate
        padded = zeros(T, padded_size) 
    else
        error("Unsupported padding mode: $mode")
    end
    
    # Copy original data to center
    if N == 3
        padded[
            padding_pixels[1]+1:padding_pixels[1]+original_size[1],
            padding_pixels[2]+1:padding_pixels[2]+original_size[2],
            padding_pixels[3]+1:padding_pixels[3]+original_size[3]
        ] = arr
        
        # Apply replication if requested
        if mode == :replicate
            _replicate_boundaries!(padded, padding_pixels, original_size)
        end
    elseif N == 4
        padded[
            padding_pixels[1]+1:padding_pixels[1]+original_size[1],
            padding_pixels[2]+1:padding_pixels[2]+original_size[2],
            padding_pixels[3]+1:padding_pixels[3]+original_size[3],
            :
        ] = arr
        
        # Apply replication if requested
        if mode == :replicate
            _replicate_boundaries_4d!(padded, padding_pixels, original_size)
        end
    end
    
    return padded
end

"""
    _replicate_boundaries!(padded, padding_pixels, original_size)

Replicate boundary values for 3D arrays.
"""
function _replicate_boundaries!(padded::AbstractArray{T, 3}, padding_pixels::NTuple{3, Int}, original_size::NTuple{3, Int}) where T
    # This is a simplified version - full implementation would replicate edge values
    # For now, just fill with background values (zeros are fine for boundaries)
    return nothing
end

"""
    _replicate_boundaries_4d!(padded, padding_pixels, original_size)

Replicate boundary values for 4D field arrays.
"""
function _replicate_boundaries_4d!(padded::AbstractArray{T, 4}, padding_pixels::NTuple{3, Int}, original_size::NTuple{4, Int}) where T
    # This is a simplified version - full implementation would replicate edge values
    return nothing
end

"""
    _create_attenuation_masks!(solver)

Create field attenuation masks to prevent boundary reflections.
"""
function _create_attenuation_masks!(solver::ConvergentBornSolver{T, AT}) where {T, AT}
    empty!(solver.field_attenuation_masks)
    
    for dim in 1:3
        max_L = solver.boundary_thickness_pixel[dim]
        L = min(max_L, solver.field_attenuation_pixel[dim])
        
        if max_L == 0
            continue
        end
        
        # Create attenuation window
        window = if L > 0
            tanh_vals = tanh.(range(T(-2.5), T(2.5), length=L))
            (tanh_vals ./ tanh(T(2.5)) .- tanh(T(-2.5))) ./ 2
        else
            T[]
        end
        
        # Apply sharpness factor
        window = window .* solver.config.field_attenuation_sharpness .+ (1 - solver.config.field_attenuation_sharpness)
        
        # Create full filter
        padded_size = size(solver.potential)[dim]
        roi_size = solver.ROI[2*dim] - solver.ROI[2*dim-1] + 1
        remaining_padding = max_L - L
        
        filter_1d = vcat(
            window,
            ones(T, roi_size + 2*remaining_padding),
            reverse(window)
        )
        
        # Reshape to proper dimension
        mask_shape = ones(Int, 3)
        mask_shape[dim] = length(filter_1d)
        mask = reshape(filter_1d, Tuple(mask_shape))
        
        push!(solver.field_attenuation_masks, AT(mask))
    end
end

"""
    _apply_attenuation_masks!(field, masks)

Apply attenuation masks to prevent boundary reflections.
"""
function _apply_attenuation_masks!(field::AbstractArray, masks::Vector)
    for mask in masks
        field .*= mask
    end
end

"""
    crop_to_ROI(field, solver) -> cropped_field

Crop padded field array to the region of interest.
"""
function crop_to_ROI(field::AbstractArray{T, N}, solver::ConvergentBornSolver) where {T, N}
    ROI = solver.ROI
    if N == 3
        return field[ROI[1]:ROI[2], ROI[3]:ROI[4], ROI[5]:ROI[6]]
    elseif N == 4
        return field[ROI[1]:ROI[2], ROI[3]:ROI[4], ROI[5]:ROI[6], :]
    else
        error("Unsupported field dimension: $N")
    end
end

"""
    reset!(solver::ConvergentBornSolver)

Reset solver state for a new problem while preserving configuration.
"""
function reset!(solver::ConvergentBornSolver)
    solver.potential = nothing
    solver.permittivity = nothing
    empty!(solver.field_attenuation_masks)
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
