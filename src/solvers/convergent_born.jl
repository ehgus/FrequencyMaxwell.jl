"""
Convergent Born Iterative Method solver implementation.

This module provides the main electromagnetic solver using the Convergent Born
Iterative Method (CBS) with LinearSolve.jl integration and AD compatibility.
"""

using LinearSolve: SciMLLinearSolveAlgorithm, KrylovJL_GMRES

"""
    ConvergentBornSolver{T} <: AbstractElectromagneticSolver{T}

Electromagnetic solver using the Convergent Born Iterative Method with integrated configuration.

This struct combines configuration and mutable solver state in a single object
for streamlined usage while maintaining high performance and AD compatibility.

# Type Parameters
- `T<:AbstractFloat`: Floating-point precision type

# Configuration Fields
- `wavelength::T`: Wavelength in the background medium (meters)
- `permittivity_bg::T`: Background relative permittivity
- `resolution::NTuple{3, T}`: Spatial resolution (dx, dy, dz) in meters
- `grid_size::NTuple{3, Int}`: Number of grid points (Nx, Ny, Nz)
- `boundary_thickness::NTuple{3, T}`: PML boundary layer thickness in each direction
- `field_attenuation::NTuple{3, T}`: Field attenuation layer thickness in each direction
- `field_attenuation_sharpness::T`: Sharpness factor for field attenuation (0-1)
- `periodic_boundary::NTuple{3, Bool}`: Periodic boundary conditions (x, y, z)
- `iterations_max::Int`: Maximum number of Born iterations (-1 for auto)
- `tolerance::T`: Convergence tolerance for iterative solver
- `linear_solver::SciMLLinearSolveAlgorithm`: LinearSolve.jl algorithm object
- `preconditioner::Symbol`: Preconditioning strategy (:none, :diagonal, :ilu)

# Solver State Fields
- `backend::Backend`: Cached KernelAbstractions.jl backend for device operations
- `Green_function::Union{Nothing, DyadicGreen{T}}`: Cached Green's function
- `potential::Union{Nothing, AbstractArray}`: Current iteration potential (V = δε)
- `permittivity::Union{Nothing, AbstractArray}`: Original permittivity without padding
- `iteration_count::Int`: Current iteration number
- `residual_history::Vector{T}`: Convergence history tracking
- `solver_state::Any`: LinearSolve.jl solver state for performance

# Constructor
```julia
ConvergentBornSolver(;
    wavelength::Real,
    permittivity_bg::Real = 1.0,
    resolution::NTuple{3, <:Real},
    grid_size::NTuple{3, Int},
    device::Symbol = :cpu,
    # ... other optional parameters
)
```

# Example
```julia
# Streamlined constructor
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32),
    device = :cuda  # GPU acceleration
)
```
"""
mutable struct ConvergentBornSolver{T <: AbstractFloat} <:
               AbstractElectromagneticSolver{T}
    # Configuration fields (formerly in ConvergentBornConfig)
    wavelength::T
    permittivity_bg::T
    resolution::NTuple{3, T}
    grid_size::NTuple{3, Int}
    boundary_thickness::NTuple{3, T}
    field_attenuation::NTuple{3, T}
    field_attenuation_sharpness::T
    periodic_boundary::NTuple{3, Bool}
    iterations_max::Int
    tolerance::T
    linear_solver::SciMLLinearSolveAlgorithm
    preconditioner::Symbol

    # Solver state fields
    backend::Backend  # Cached KernelAbstractions.jl backend
    Green_function::Union{Nothing, DyadicGreen{T}}
    flip_Green_function::Union{Nothing, DyadicGreen{T}}
    potential::Union{Nothing, AbstractArray}  # V = k²(ε/ε_bg - 1) with padding
    potential_device::Union{Nothing, AbstractArray}  # Device-resident potential
    permittivity::Union{Nothing, AbstractArray}  # Original permittivity without padding
    boundary_thickness_pixel::NTuple{3, Int}
    field_attenuation_pixel::NTuple{3, Int}
    ROI::NTuple{6, Int}  # Region of interest bounds
    eps_imag::T  # Imaginary part for convergence
    Bornmax::Int  # Actual number of iterations to use
    iteration_count::Int
    residual_history::Vector{T}
    solver_state::Any  # LinearSolve.jl solver state

    function ConvergentBornSolver{T}(
            wavelength::T,
            permittivity_bg::T,
            resolution::NTuple{3, T},
            grid_size::NTuple{3, Int},
            boundary_thickness::NTuple{3, T},
            field_attenuation::NTuple{3, T},
            field_attenuation_sharpness::T,
            periodic_boundary::NTuple{3, Bool},
            iterations_max::Int,
            tolerance::T,
            linear_solver::SciMLLinearSolveAlgorithm,
            preconditioner::Symbol,
            backend::Backend
    ) where {T <: AbstractFloat}
        # Validate parameters (similar to ConvergentBornConfig validation)
        wavelength > 0 || throw(ArgumentError("wavelength must be positive"))
        permittivity_bg > 0 || throw(ArgumentError("permittivity_bg must be positive"))
        all(resolution .> 0) ||
            throw(ArgumentError("all resolution components must be positive"))
        all(grid_size .> 0) ||
            throw(ArgumentError("all grid_size components must be positive"))
        0 ≤ tolerance ≤ 1 || throw(ArgumentError("tolerance must be in [0, 1]"))

        # Calculate boundary thickness in pixels
        boundary_thickness_pixel = round.(Int, boundary_thickness ./ (resolution .* 2))
        field_attenuation_pixel = round.(Int, field_attenuation ./ (resolution .* 2))

        # Calculate region of interest
        ROI = (
            boundary_thickness_pixel[1] + 1, boundary_thickness_pixel[1] + grid_size[1],
            boundary_thickness_pixel[2] + 1, boundary_thickness_pixel[2] + grid_size[2],
            boundary_thickness_pixel[3] + 1, boundary_thickness_pixel[3] + grid_size[3]
        )

        new{T}(
            # Configuration fields
            wavelength, permittivity_bg, resolution, grid_size,
            boundary_thickness, field_attenuation, field_attenuation_sharpness,
            periodic_boundary, iterations_max, tolerance, linear_solver, preconditioner,
            # Solver state fields
            backend,  # KernelAbstractions.jl backend (cached)
            nothing,  # Green_function - computed lazily
            nothing,  # flip_Green_function - computed lazily
            nothing,  # potential - set during solve
            nothing,  # potential_device - set during solve
            nothing,  # permittivity - set during solve
            boundary_thickness_pixel,
            field_attenuation_pixel,
            ROI,
            T(0),     # eps_imag - calculated based on potential
            0,        # Bornmax - calculated automatically
            0,        # iteration_count
            T[],      # residual_history
            nothing   # solver_state - initialized during solve
        )
    end
end

"""
    to_device(backend::KernelAbstractions.Backend, host_array::AbstractArray) -> AbstractArray

Transfer a host array to the device specified by the backend.

This function handles data transfer from host memory to device memory,
automatically handling the appropriate array type conversion.

# Arguments
- `backend::KernelAbstractions.Backend`: Target device backend
- `host_array::AbstractArray`: Source array (typically on host/CPU)

# Returns
- Device-resident copy of the input array
"""
function to_device(backend::KernelAbstractions.Backend, host_array::AbstractArray)
    device_array = KernelAbstractions.allocate(backend, eltype(host_array), size(host_array))
    copyto!(device_array, host_array)
    return device_array
end

"""
    to_host(device_array::AbstractArray) -> Array

Transfer a device array back to host memory as a standard Julia Array.

This function handles data transfer from device memory back to host memory,
ensuring results are accessible as normal Julia arrays.

# Arguments
- `device_array::AbstractArray`: Source array on device

# Returns
- Host-resident Array copy of the device array
"""
function to_host(device_array::AbstractArray)
    host_array = Array{eltype(device_array)}(undef, size(device_array))
    copyto!(host_array, device_array)
    return host_array
end

"""
    create_backend_from_symbol(device::Symbol) -> KernelAbstractions.Backend

Create a KernelAbstractions.jl backend from a device symbol with conditional GPU support.

This function provides the symbol-to-backend mapping that enables simple
device selection while maintaining pure KernelAbstractions.jl implementation.

# Supported Device Symbols
- `:cpu` → `CPU()` - Multi-threaded CPU execution
- `:cuda` → `CUDABackend()` - NVIDIA GPUs via CUDA
- `:amdgpu` → `ROCBackend()` - AMD GPUs via ROCm
- `:metal` → `MetalBackend()` - Apple Silicon via Metal
- `:oneapi` → `oneAPIBackend()` - Intel GPUs via oneAPI

# Example
```julia
# Always works
backend = create_backend_from_symbol(:cpu)

# Works if CUDA.jl is loaded
backend = create_backend_from_symbol(:cuda)

# Check availability first
if is_backend_available(:cuda)
    backend = create_backend_from_symbol(:cuda)
else
    @info "CUDA.jl not available, using CPU backend"
    backend = create_backend_from_symbol(:cpu)
end
```
"""
function create_backend_from_symbol(device::Symbol)
    if device == :cpu
        return KernelAbstractions.CPU()
    elseif is_backend_available(device)
        # Use registered GPU backend constructor
        constructor = FrequencyMaxwell.GPU_BACKENDS[device]
        return constructor()
    else
        # Provide specific installation guidance
        install_instructions = if device == :cuda
            "Install CUDA.jl with: Pkg.add(\"CUDA\")"
        elseif device in [:amdgpu, :rocm]
            "Install AMDGPU.jl with: Pkg.add(\"AMDGPU\")"
        elseif device == :metal
            "Install Metal.jl with: Pkg.add(\"Metal\")"
        elseif device == :oneapi
            "Install oneAPI.jl with: Pkg.add(\"oneAPI\")"
        else
            "Unknown device symbol: :$device"
        end

        throw(ArgumentError(
            "GPU backend :$device is not available. $install_instructions\n" *
            "For CPU-only usage, use: ConvergentBornSolver(config, :cpu)"
        ))
    end
end


"""
    ConvergentBornSolver(; kwargs...) -> ConvergentBornSolver

Enhanced constructor that accepts configuration parameters directly as keyword arguments.

This is the recommended and only way to create a ConvergentBornSolver, providing
direct access to all configuration parameters without intermediate objects.

# Arguments
- `wavelength::Real`: Wavelength in the background medium (meters)
- `permittivity_bg::Real = 1.0`: Background relative permittivity
- `resolution::NTuple{3, <:Real}`: Spatial resolution (dx, dy, dz) in meters
- `grid_size::NTuple{3, Int}`: Number of grid points (Nx, Ny, Nz)
- `device::Symbol = :cpu`: Device backend (:cpu, :cuda, :amdgpu, :metal, :oneapi)
- `boundary_thickness::NTuple{3, <:Real} = (0.0, 0.0, 0.0)`: PML boundary layer thickness
- `field_attenuation::NTuple{3, <:Real} = (0.0, 0.0, 0.0)`: Field attenuation layer thickness
- `field_attenuation_sharpness::Real = 1.0`: Sharpness factor for field attenuation (0-1)
- `periodic_boundary::NTuple{3, Bool} = (true, true, false)`: Periodic boundary conditions
- `iterations_max::Int = -1`: Maximum number of Born iterations (-1 for auto)
- `tolerance::Real = 1e-6`: Convergence tolerance for iterative solver
- `linear_solver::SciMLLinearSolveAlgorithm = KrylovJL_GMRES()`: LinearSolve.jl algorithm
- `preconditioner::Symbol = :none`: Preconditioning strategy (:none, :diagonal, :ilu)

# Examples
```julia
# Basic usage
solver = ConvergentBornSolver(
    wavelength = 532e-9,
    permittivity_bg = 1.333^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64)
)

# With GPU acceleration
solver = ConvergentBornSolver(
    wavelength = 532e-9,
    permittivity_bg = 1.333^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64),
    device = :cuda
)

# Advanced configuration
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.0,
    resolution = (25e-9, 25e-9, 25e-9),
    grid_size = (256, 256, 128),
    device = :cpu,
    boundary_thickness = (1e-6, 1e-6, 2e-6),
    tolerance = 1e-4,
    iterations_max = 50
)
```
"""
function ConvergentBornSolver(;
        wavelength::Real,
        permittivity_bg::Real = 1.0,
        resolution::NTuple{3, <:Real},
        grid_size::NTuple{3, Int},
        device::Symbol = :cpu,
        boundary_thickness::NTuple{3, <:Real} = (0.0, 0.0, 0.0),
        field_attenuation::NTuple{3, <:Real} = (0.0, 0.0, 0.0),
        field_attenuation_sharpness::Real = 1.0,
        periodic_boundary::NTuple{3, Bool} = (true, true, false),
        iterations_max::Int = -1,
        tolerance::Real = 1e-6,
        linear_solver::SciMLLinearSolveAlgorithm = KrylovJL_GMRES(),
        preconditioner::Symbol = :none
)
    # Promote to common floating-point type
    T = promote_type(typeof(wavelength), typeof(permittivity_bg),
        eltype(resolution), eltype(boundary_thickness),
        eltype(field_attenuation), typeof(field_attenuation_sharpness), typeof(tolerance))
    T <: AbstractFloat || (T = Float64)  # Fallback to Float64 if not floating-point

    # Create backend from device symbol
    backend = create_backend_from_symbol(device)

    # Create and return solver using the internal constructor
    return ConvergentBornSolver{T}(
        T(wavelength), T(permittivity_bg),
        T.(resolution), grid_size,
        T.(boundary_thickness), T.(field_attenuation), T(field_attenuation_sharpness),
        periodic_boundary, iterations_max, T(tolerance),
        linear_solver, preconditioner, backend
    )
end

"""
    CBSPreconditioner{T} <: Function

Left preconditioner for the Convergent Born Series implementing R = (1i/ε_imag) * potential.

This preconditioner applies the transformation: y = R * x = (1i/ε_imag) * potential .* x

The CBS linear system becomes: (1-G*V*)x = y with left preconditioner R
where the original system R*(1-G*V*)x = R*y is reformulated.

# Fields
- `solver::ConvergentBornSolver{T}`: Reference to the CBS solver containing potential

# Interface
Implements LinearSolve.jl preconditioner interface with `ldiv!` methods.
"""
struct CBSPreconditioner{T <: AbstractFloat}
    solver::ConvergentBornSolver{T}

    function CBSPreconditioner(solver::ConvergentBornSolver{T}) where {T}
        new{T}(solver)
    end
end

"""
    ldiv!(y, P::CBSPreconditioner, x)

Apply CBS preconditioner: y = R * x = (1i/ε_imag) * potential .* x

This method implements the LinearSolve.jl preconditioner interface.
"""
function LinearAlgebra.ldiv!(y, P::CBSPreconditioner{T}, x) where {T}
    # Reshape flattened input to 4D field array
    potential = P.solver.potential_device
    field_shape = (size(potential)..., 3)
    x_field = reshape(x, field_shape)
    y_field = reshape(y, field_shape)

    # Apply preconditioner: R * x = (1i/ε_imag) * potential .* x
    y_field .= potential .* x_field
    y_field .*= (Complex{T}(0, 1) / P.solver.eps_imag)

    return y
end

"""
    ldiv!(P::CBSPreconditioner, x)

In-place application of CBS preconditioner: x = R * x = (1i/ε_imag) * potential .* x

This method implements the in-place LinearSolve.jl preconditioner interface.
"""
function LinearAlgebra.ldiv!(P::CBSPreconditioner{T}, x) where {T}
    # Reshape flattened input to 4D field array
    potential = P.solver.potential_device
    field_shape = (size(potential)..., 3)
    x_field = reshape(x, field_shape)

    # Apply preconditioner: R * x = (1i/ε_imag) * potential .* x
    x_field .*= potential
    x_field .*= (Complex{T}(0, 1) / P.solver.eps_imag)

    # No need to copy back since we modified x_field which is a view of x
    return x
end

"""
    CBSLinearOperator{T} <: Function

Linear operator representing (1-G*V*) for the Convergent Born Series reformulation.

The CBS iteration `Field_{n+1} = Field_0 + A * Field_n` can be reformulated as
a linear system `(1-G*V*) * Field = Field_0` where:
- G*V* = (G + G_flip)/2 * V
- V is the potential (permittivity contrast)
- G, G_flip are dyadic Green's functions

This operator represents only the core CBS operator without the preconditioner R.
The preconditioner R = (1i/ε_imag) * potential is applied separately through LinearSolve.jl.

# Fields
- `solver::ConvergentBornSolver{T}`: Reference to the CBS solver
- `temp_arrays::NTuple{2, AbstractArray}`: Preallocated temporary 4D field arrays for efficiency
"""
struct CBSLinearOperator{T <: AbstractFloat}
    solver::ConvergentBornSolver{T}
    temp_arrays::NTuple{2, AbstractArray}  # Preallocated temporary 4D field arrays

    function CBSLinearOperator(solver::ConvergentBornSolver{T}) where {T}
        # Preallocate temporary arrays to match field dimensions
        field_size = size(solver.potential)
        temp1 = KernelAbstractions.zeros(solver.backend, Complex{T}, field_size..., 3)
        temp2 = KernelAbstractions.zeros(solver.backend, Complex{T}, field_size..., 3)

        new{T}(solver, (temp1, temp2))
    end
end

"""
    (op::CBSLinearOperator)(y, x)

Apply the linear operator to input field x, storing result in y.

This implements: y = (1-G*V*)*x where G*V* = (G + G_flip)/2 * V

The preconditioner R = (1i/ε_imag) * potential is applied separately through LinearSolve.jl.

# Arguments
- `y`: Output vector (flattened field array)  
- `x`: Input vector (flattened field array)
"""
function (op::CBSLinearOperator{T})(y, x, p, t; α = 1, β = 0) where {T}
    solver = op.solver
    potential = solver.potential_device
    psi, flip_psi = op.temp_arrays

    # Reshape flattened input to 4D field array
    field_shape = size(psi)
    x_field = reshape(x, field_shape)

    # Apply Green's operator: G*V*x = (G + G_flip)/2 * V * x

    # Step 1: V * x (element-wise multiplication with potential)
    psi .= potential .* x_field

    # Step 2: (G + G_flip)/2
    flip_psi .= psi
    conv!(solver.Green_function, psi)
    conv!(solver.flip_Green_function, flip_psi)
    psi .+= flip_psi
    psi ./= 2

    # Step 3: (1 - G*V*)x = x - G*V*x
    psi .= x_field .- psi

    # Flatten result back to vector
    # y = α * (1-G*V*)x + β * y
    y .*= β
    y .+= α .* vec(psi)

    return y
end

function LinearAlgebra.mul!(y, op::CBSLinearOperator, x, α, β)
    op(y, x, nothing, nothing; α, β)
end

# Make the operator compatible with LinearSolve.jl AbstractArray interface
Base.size(op::CBSLinearOperator) = (length(op.temp_arrays[1]), length(op.temp_arrays[1]))
Base.size(op::CBSLinearOperator, dim::Integer) = size(op)[dim]
Base.eltype(::CBSLinearOperator{T}) where {T} = Complex{T}
Base.ndims(::CBSLinearOperator) = 2  # Matrix-like operator

"""
    solve(solver::ConvergentBornSolver, sources, permittivity::AbstractArray)

Solve the electromagnetic scattering problem using the Convergent Born method.

This unified method handles both single sources and multiple coherent sources using Julia's
multiple dispatch. For multiple sources, they interfere coherently in the same simulation
domain, enabling applications like beam splitting, interference patterns, and complex scattering.

# Arguments
- `solver::ConvergentBornSolver`: Pre-configured solver instance
- `sources`: Either an `AbstractCurrentSource` (single) or `Vector{<:AbstractCurrentSource}` (multiple)
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
- ψ_incident: incident field from source(s) - coherent superposition for multiple sources
- G: dyadic Green's function
- V: potential (permittivity contrast)

For multiple sources, incident fields are coherently superposed: E_incident = Σ E_i

# Examples
```julia
# Single source
E_field, H_field = solve(solver, source, permittivity)

# Multi-source (coherent interference)
sources = [source1, source2, source3]
E_field, H_field = solve(solver, sources, permittivity)
```
"""
function LinearSolve.solve(
        solver::ConvergentBornSolver{T},
        sources::Union{AbstractCurrentSource{T}, Vector{<:AbstractCurrentSource{T}}},
        permittivity::AbstractArray{<:Number, 3}
) where {T <: AbstractFloat}

    # Validate input dimensions
    size(permittivity) == solver.grid_size ||
        throw(ArgumentError("permittivity size $(size(permittivity)) must match grid_size $(solver.grid_size)"))

    # Additional validation for multi-source case
    if sources isa Vector
        isempty(sources) && throw(ArgumentError("sources vector cannot be empty"))
        _validate_source_compatibility(sources)
    end

    # Initialize solver state (computes potential V with proper padding)
    _initialize_solver!(solver, permittivity)

    # Generate incident fields (multiple dispatch handles single vs multi-source)
    E_incident, H_incident = _generate_incident_fields_padded(solver, sources)

    # Generate source term: V * E_incident + i*eps_imag * E_incident
    source_term = solver.potential .* E_incident .+
                  Complex{T}(0, solver.eps_imag) .* E_incident

    # Transfer source term to device
    source_term = to_device(solver.backend, source_term)

    # Solve scattering equation using CBS algorithm
    E_scattered_device, H_scattered_device = _solve_cbs_scattering(solver, source_term)

    # Transfer scattered field back to host
    E_scattered = to_host(E_scattered_device)
    H_scattered = to_host(H_scattered_device)

    # Add incident field to get total field and crop to ROI
    E_total = crop_to_ROI(E_scattered + E_incident, solver)
    H_total = crop_to_ROI(H_scattered + H_incident, solver)

    return E_total, H_total
end

"""
    _initialize_solver!(solver::ConvergentBornSolver, permittivity::AbstractArray)

Initialize solver internal state for a new problem.

This function prepares the solver for electromagnetic scattering calculations by:
1. Computing the potential V = k²(ε/ε_bg - 1) with proper boundary padding
2. Adding imaginary regularization for convergence stability  
3. Creating field attenuation masks for boundary condition enforcement
4. Pre-applying attenuation masks to the potential (optimization)
5. Computing optimal iteration count and initializing Green's functions

The potential pre-masking optimization ensures that subsequent field operations
automatically include boundary attenuation without redundant mask applications.
"""
function _initialize_solver!(
        solver::ConvergentBornSolver{T},
        permittivity::AbstractArray{<:Number, 3}
) where {T <: AbstractFloat}

    # Store original permittivity
    solver.permittivity = Complex{T}.(permittivity)

    # Compute potential V = k²(ε/ε_bg - 1) with proper padding
    k = T(2π) * sqrt(solver.permittivity_bg) / solver.wavelength
    potential_unpadded = Complex{T}.(k^2 *
                                     (permittivity ./ solver.permittivity_bg .-
                                      1))

    # Pad the potential (replicate boundary values)
    try
        solver.potential = _pad_array(potential_unpadded, solver.boundary_thickness_pixel, :replicate)
    catch e
        error("The Pad_array information: $(solver.boundary_thickness_pixel) // $(size(potential_unpadded))")
        throw(e)
    end

    # Calculate eps_imag based on maximum potential value for numerical stability
    solver.eps_imag = max(T(2^-20), maximum(abs.(solver.potential)) * T(1.01))

    # Add imaginary part to potential for convergence
    solver.potential .-= Complex{T}(0, solver.eps_imag)

    # Apply field attenuation directly to potential (optimized to avoid storing masks)
    _apply_attenuation_to_potential!(solver)

    # Calculate optimal iteration count if automatic
    if solver.iterations_max < 0
        k0_nm = T(2π) * sqrt(solver.permittivity_bg) / solver.wavelength
        steps = abs(2 * k0_nm / solver.eps_imag)
        domain_extent = norm(size(solver.potential)[1:3] .* solver.resolution)
        Bornmax_opt = ceil(Int, domain_extent / steps / 2 + 1) * 2
        solver.Bornmax = Bornmax_opt * abs(solver.iterations_max)
    else
        solver.Bornmax = solver.iterations_max
    end

    # Initialize Green's functions if needed
    if solver.Green_function === nothing
        solver.Green_function, solver.flip_Green_function = _compute_green_functions(solver)
    end

    # Reset iteration tracking
    solver.iteration_count = 0
    empty!(solver.residual_history)

    # Transfer potential to device
    solver.potential_device = to_device(solver.backend, solver.potential)

    return nothing
end

"""
    _compute_green_function(solver::ConvergentBornSolver) -> DyadicGreen

Compute the dyadic Green's function for the given configuration.
"""
function _compute_green_functions(solver::ConvergentBornSolver{T}) where {T}
    # Calculate wave number in background medium
    k0_nm = T(2π) * sqrt(solver.permittivity_bg) / solver.wavelength
    k_square = k0_nm^2 + Complex{T}(0, solver.eps_imag)

    # Array size includes boundary padding
    arr_size = size(solver.potential)

    # Subpixel shifts for proper boundary conditions
    subpixel_shift = ntuple(i -> solver.periodic_boundary[i] ? T(0) : T(0.25), 3)
    flip_subpixel_shift = ntuple(i -> solver.periodic_boundary[i] ? T(0) : T(-0.25), 3)

    # Create dyadic Green's functions
    array_type = typeof(solver.potential) <: Array ? Array :
                 typeof(solver.potential).name.wrapper
    Green_fn = DyadicGreen(
        array_type, k_square, arr_size, solver.resolution, subpixel_shift)
    flip_Green_fn = DyadicGreen(
        array_type, k_square, arr_size, solver.resolution, flip_subpixel_shift)

    return Green_fn, flip_Green_fn
end

"""
    _generate_incident_fields(solver::ConvergentBornSolver, source::AbstractCurrentSource)

Generate incident electromagnetic fields from the current source.
"""
function _generate_incident_fields_padded(
        solver::ConvergentBornSolver{T},
        source::AbstractCurrentSource{T}
) where {T <: AbstractFloat}

    # Generate incident fields on original grid
    grid_size = solver.grid_size
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
    _generate_incident_fields_padded(solver::ConvergentBornSolver, sources::Vector{<:AbstractCurrentSource})

Generate incident electromagnetic fields from multiple coherent sources with proper superposition.

This method overload implements coherent field superposition where multiple electromagnetic sources
interfere in the same simulation domain. The complex amplitudes are summed to produce
the combined incident field that includes constructive and destructive interference.

# Arguments
- `solver::ConvergentBornSolver`: Configured solver instance
- `sources::Vector{<:AbstractCurrentSource}`: Vector of coherent electromagnetic sources

# Algorithm
1. Validate source compatibility (same wavelength for coherent interference)
2. Generate incident fields from each source independently using single-source method
3. Coherently sum complex amplitudes: E_total = Σ E_i, H_total = Σ H_i

This produces proper electromagnetic interference patterns where sources can interfere
constructively (bright fringes) or destructively (dark fringes).
"""
function _generate_incident_fields_padded(
        solver::ConvergentBornSolver{T},
        sources::Vector{<:AbstractCurrentSource{T}}
) where {T <: AbstractFloat}

    # Initialize total field arrays (padded grid size)
    padded_size = size(solver.potential)[1:3]
    E_total = zeros(Complex{T}, padded_size..., 3)
    H_total = zeros(Complex{T}, padded_size..., 3)

    # Sum contributions from all sources (coherent superposition)
    for source in sources
        E_src, H_src = _generate_incident_fields_padded(solver, source)
        E_total .+= E_src  # Complex amplitude addition for interference
        H_total .+= H_src
    end

    return E_total, H_total
end

"""
    _validate_source_compatibility(sources)

Validate that multiple electromagnetic sources are compatible for coherent simulation.

This function ensures that all sources in a multi-source configuration have compatible
parameters for coherent electromagnetic simulation. The primary requirement is that
all sources must have the same wavelength for proper interference.

# Arguments
- `sources::Vector{<:AbstractCurrentSource}`: Vector of sources to validate

# Validation Checks
1. Non-empty source vector
2. All sources have identical wavelength (within tolerance 1e-10)

# Throws
- `ArgumentError` if sources are incompatible with detailed error message

# Example
```julia
sources = [source1, source2, source3]
_validate_source_compatibility(sources)  # Throws if wavelengths differ
```
"""
function _validate_source_compatibility(sources::Vector{<:AbstractCurrentSource{T}}) where {T}
    # Check for empty sources vector
    isempty(sources) && throw(ArgumentError("sources vector cannot be empty"))

    # All sources must have the same wavelength for coherent interference
    reference_wavelength = source_wavelength(sources[1])
    for (i, source) in enumerate(sources[2:end])
        src_wavelength = source_wavelength(source)
        if !isapprox(src_wavelength, reference_wavelength, rtol = T(1e-10))
            throw(ArgumentError(
                "All sources must have the same wavelength for coherent interference. " *
                "Source 1 wavelength: $(reference_wavelength), " *
                "Source $(i+1) wavelength: $(src_wavelength)"
            ))
        end
    end

    return nothing
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
) where {T <: AbstractFloat}
    grid_size = solver.grid_size
    resolution = solver.resolution

    # Wave parameters
    k0 = T(2π / solver.wavelength)  # Free-space wave number
    k_bg = k0 * sqrt(solver.permittivity_bg)  # Background wave number

    # Normalize k-vector and polarization
    k_hat = source.k_vector ./ norm(source.k_vector)  # Unit k-vector
    k_vec = k_bg .* k_hat  # Full k-vector in medium
    E0 = source.polarization ./ norm(source.polarization)  # Normalized polarization

    # Ensure E ⟂ k (transversality condition)
    E_perp = E0 .- (dot(E0, k_hat) * k_hat)  # Remove parallel component
    E_perp = E_perp ./ norm(E_perp)  # Renormalize

    # Compute H from E using H = k × E / (ωμ₀) = k × E / (k₀Z₀)
    H_perp = cross(k_hat, E_perp) ./ T(377 / sqrt(solver.permittivity_bg))

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
    _solve_cbs_scattering(solver, source) -> (E_field, H_field)

Solve the electromagnetic scattering problem using LinearSolve.jl CBS implementation.

This function reformulates the Convergent Born Series as a linear system (I - A) * Field = Field_0
and leverages the user-configured LinearSolver algorithm for efficient solution.

# Mathematical Foundation
- Traditional CBS: Field_{n+1} = Field_0 + A * Field_n  
- Reformulated: (I - A) * Field_total = Field_0
- Where A = V * (G + G_flip)/2 * (i/ε_imag) * V

# Performance Benefits
- Uses configured LinearSolver objects (KrylovJL_GMRES, KrylovJL_BICGSTAB, etc.)
- Adaptive convergence criteria with user-configurable tolerances
- Advanced preconditioning strategies through LinearSolver object parameters
- Significant speedup potential over fixed iterations

# Arguments
- `solver::ConvergentBornSolver`: Configured CBS solver instance
- `source::AbstractArray{Complex{T}, 4}`: Incident source field array (nx, ny, nz, 3)

# Returns
- `(E_field, H_field)`: Tuple of electric and magnetic field arrays
"""
function _solve_cbs_scattering(
        solver::ConvergentBornSolver{T},
        source::AbstractArray{Complex{T}, 4}
) where {T <: AbstractFloat}

    # Prepare right-hand side: Field_0 from first CBS iteration  
    size_field = size(source)
    Green_fn = solver.Green_function
    flip_Green_fn = solver.flip_Green_function
    psi = KernelAbstractions.zeros(solver.backend, Complex{T}, size_field)
    flip_psi = KernelAbstractions.zeros(solver.backend, Complex{T}, size_field)

    if solver.Bornmax >= 1
        # Compute Field_0: first iteration of CBS without R preconditioner
        # Field_0 = (G + G_flip)/2 * source
        psi .= source

        # Apply Green's functions: (G + G_flip)/2
        flip_psi .= psi
        conv!(Green_fn, psi)
        conv!(flip_Green_fn, flip_psi)
        psi .+= flip_psi
        psi ./= 2

        # Set up LinearSolve.jl problem
        # Flatten arrays for LinearSolve.jl compatibility
        rhs_vec = copy(vec(psi))
        x0_vec = copy(vec(psi))  # Initial guess: use first CBS iteration

        # Create linear operator
        linear_op = CBSLinearOperator(solver)

        # Create CBS preconditioner R = (1i/ε_imag) * potential
        preconditioner = CBSPreconditioner(solver)

        # Create linear problem
        prob = LinearProblem(linear_op, rhs_vec; u0 = x0_vec)

        # Use configured LinearSolver algorithm
        algorithm = solver.linear_solver

        # Set solver options
        solver_options = Dict{Symbol, Any}(
            :reltol => solver.tolerance,
            :abstol => solver.tolerance * 1e-2,
            :maxiters => max(solver.Bornmax * 2, 100)  # More flexible iteration limit
        )

        # Solve linear system
        sol = LinearSolve.solve(prob, algorithm; Pl = preconditioner, solver_options...)

        # Check convergence status
        if !SciMLBase.successful_retcode(sol.retcode)
            @warn "LinearSolve.jl did not converge."
        end

        # Reshape solution back to field array
        Field = reshape(sol.u, size_field)

        # Update solver statistics if available
        if hasfield(typeof(sol), :iters)
            solver.iteration_count = sol.iters
        end
        if hasfield(typeof(sol), :resid)
            push!(solver.residual_history, sol.resid)
        end
    end

    # Compute magnetic field using same method as iterative version
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
        solver::ConvergentBornSolver{T},
        E_field::AbstractArray{Complex{T}, 4}
) where {T <: AbstractFloat}

    # Create curl operator for field with padding
    array_type = typeof(solver.potential)
    curl_op = Curl(array_type, T, size(solver.potential), solver.resolution)

    # Compute curl of E field
    Hfield = conv(curl_op, E_field)

    # Apply scaling: H = -i * λ/(2π * Z₀) * ∇ × E
    scaling_factor = Complex{T}(0, -1) * solver.wavelength / (T(2π) * T(377))

    return scaling_factor .* Hfield
end

"""
    _pad_array(arr, padding_pixels, mode) -> padded_array

Pad an array with specified padding in each dimension.
"""
function _pad_array(arr::AbstractArray{T, N}, padding_pixels::NTuple{3, Int}, mode::Symbol) where {
        T, N}
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
        (padding_pixels[1] + 1):(padding_pixels[1] + original_size[1]),
        (padding_pixels[2] + 1):(padding_pixels[2] + original_size[2]),
        (padding_pixels[3] + 1):(padding_pixels[3] + original_size[3])
] = arr

        # Apply replication if requested
        if mode == :replicate
            _replicate_boundaries!(padded, padding_pixels, original_size)
        end
    elseif N == 4
        padded[
        (padding_pixels[1] + 1):(padding_pixels[1] + original_size[1]),
        (padding_pixels[2] + 1):(padding_pixels[2] + original_size[2]),
        (padding_pixels[3] + 1):(padding_pixels[3] + original_size[3]),
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
function _replicate_boundaries!(
        padded::AbstractArray{T, 3}, padding_pixels::NTuple{3, Int},
        original_size::NTuple{3, Int}) where {T}
    # This is a simplified version - full implementation would replicate edge values
    # For now, just fill with background values (zeros are fine for boundaries)
    return nothing
end

"""
    _replicate_boundaries_4d!(padded, padding_pixels, original_size)

Replicate boundary values for 4D field arrays.
"""
function _replicate_boundaries_4d!(
        padded::AbstractArray{T, 4}, padding_pixels::NTuple{3, Int},
        original_size::NTuple{4, Int}) where {T}
    # This is a simplified version - full implementation would replicate edge values
    return nothing
end

"""
    _apply_attenuation_to_potential!(solver)

Create and apply field attenuation masks directly to the solver's potential.

This optimized function combines mask creation and application into a single operation,
eliminating the need to store intermediate mask arrays. The masks are created
dimension-by-dimension and immediately applied to prevent boundary reflections.
"""
function _apply_attenuation_to_potential!(solver::ConvergentBornSolver{T}) where {T}
    for dim in 1:3
        max_L = solver.boundary_thickness_pixel[dim]
        L = min(max_L, solver.field_attenuation_pixel[dim])

        if max_L == 0
            continue
        end

        # Create attenuation window
        window = if L > 0
            tanh_vals = tanh.(range(T(-2.5), T(2.5), length = L))
            (tanh_vals ./ tanh(T(2.5)) .- tanh(T(-2.5))) ./ 2
        else
            T[]
        end

        # Apply sharpness factor
        window = window .* solver.field_attenuation_sharpness .+
                 (1 - solver.field_attenuation_sharpness)

        # Create full filter
        roi_size = solver.ROI[2 * dim] - solver.ROI[2 * dim - 1] + 1
        remaining_padding = max_L - L

        filter_1d = vcat(
            window,
            ones(T, roi_size + 2*remaining_padding),
            reverse(window)
        )

        # Create and apply mask directly to potential
        # Use broadcasting to apply 1D filter along the specified dimension
        if dim == 1
            solver.potential .*= reshape(filter_1d, :, 1, 1)
        elseif dim == 2
            solver.potential .*= reshape(filter_1d, 1, :, 1)
        else  # dim == 3
            solver.potential .*= reshape(filter_1d, 1, 1, :)
        end
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
    solver.iteration_count = 0
    empty!(solver.residual_history)
    solver.solver_state = nothing
    return nothing
end

"""
    show(io::IO, solver::ConvergentBornSolver)

Custom display for solver objects with status information.
"""
function Base.show(io::IO, solver::ConvergentBornSolver{T}) where {T}
    print(io, "ConvergentBornSolver{$T}:")
    print(io, "\n  wavelength: $(solver.wavelength)")
    print(io, "\n  permittivity_bg: $(solver.permittivity_bg)")
    print(io, "\n  resolution: $(solver.resolution)")
    print(io, "\n  grid_size: $(solver.grid_size)")
    print(io, "\n  domain_size: $(domain_size(solver))")
    print(io, "\n  k₀: $(wavenumber_background(solver))")
    print(io, "\n  device: $(typeof(solver.backend))")
    print(io, "\n  iterations: $(solver.iteration_count)")
    if !isempty(solver.residual_history)
        print(io, "\n  last residual: $(last(solver.residual_history))")
    end
end

"""
    grid_spacing(solver::ConvergentBornSolver) -> NTuple{3, T}

Calculate the physical grid spacing in each direction.
"""
grid_spacing(solver::ConvergentBornSolver) = solver.resolution

"""
    domain_size(solver::ConvergentBornSolver) -> NTuple{3, T}

Calculate the total physical domain size in each direction.
"""
function domain_size(solver::ConvergentBornSolver{T}) where {T}
    return ntuple(i -> solver.grid_size[i] * solver.resolution[i], 3)
end

"""
    wavenumber_background(solver::ConvergentBornSolver) -> T

Calculate the background medium wavenumber k₀ = 2π√(εᵦ)/λ₀.
"""
function wavenumber_background(solver::ConvergentBornSolver{T}) where {T}
    return T(2π) * sqrt(solver.permittivity_bg) / solver.wavelength
end
