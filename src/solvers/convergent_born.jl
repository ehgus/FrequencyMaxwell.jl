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
- `boundary_conditions::NTuple{3, AbstractBoundaryCondition{T}}`: Boundary conditions per dimension (x, y, z)
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
    boundary_conditions::Union{AbstractBoundaryCondition, NTuple{3, <:AbstractBoundaryCondition}},
    device::Symbol = :cpu,
    # ... other optional parameters
)
```

# Examples
```julia
# New interface with explicit boundary conditions
bc_absorbing = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6,
    sharpness = 0.9
)
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32),
    boundary_conditions = bc_absorbing,  # Same for all dimensions
    device = :cuda
)

# Mixed boundary conditions (periodic in x,y, absorbing in z)
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    boundary_conditions = (
        PeriodicBoundaryCondition(),
        PeriodicBoundaryCondition(),
        bc_absorbing
    )
)
```

# Boundary Condition Integration

The solver now uses the boundary condition abstraction:
- Padding pixels calculated from boundary conditions via `padding_pixels()`
- Subpixel shifts extracted via `subpixel_shift()`
- Field attenuation applied via `apply_attenuation!()`
- Green's function averaging controlled via `requires_averaging()`

This eliminates hardcoded boundary logic and enables flexible boundary configurations.
"""
mutable struct ConvergentBornSolver{T <: AbstractFloat} <:
               AbstractElectromagneticSolver{T}
    # Configuration fields
    wavelength::T
    permittivity_bg::T
    resolution::NTuple{3, T}
    grid_size::NTuple{3, Int}
    boundary_conditions::NTuple{3, AbstractBoundaryCondition{T}}
    iterations_max::Int
    tolerance::T
    linear_solver::SciMLLinearSolveAlgorithm
    preconditioner::Symbol

    # Solver state fields
    backend::Backend  # Cached KernelAbstractions.jl backend
    Green_function::Union{Nothing, DyadicGreen{T}}
    potential::Union{Nothing, AbstractArray}  # V = k²(ε/ε_bg - 1) with padding
    potential_device::Union{Nothing, AbstractArray}  # Device-resident potential
    permittivity::Union{Nothing, AbstractArray}  # Original permittivity without padding
    boundary_thickness_pixel::NTuple{3, Int}  # Derived from boundary conditions
    field_attenuation_pixel::NTuple{3, Int}   # Derived from boundary conditions
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
            boundary_conditions::NTuple{3, AbstractBoundaryCondition{T}},
            iterations_max::Int,
            tolerance::T,
            linear_solver::SciMLLinearSolveAlgorithm,
            preconditioner::Symbol,
            backend::Backend
    ) where {T <: AbstractFloat}
        # Validate parameters
        wavelength > 0 || throw(ArgumentError("wavelength must be positive"))
        permittivity_bg > 0 || throw(ArgumentError("permittivity_bg must be positive"))
        all(resolution .> 0) ||
            throw(ArgumentError("all resolution components must be positive"))
        all(grid_size .> 0) ||
            throw(ArgumentError("all grid_size components must be positive"))
        0 ≤ tolerance ≤ 1 || throw(ArgumentError("tolerance must be in [0, 1]"))

        # Calculate boundary thickness in pixels from boundary conditions
        boundary_thickness_pixel = ntuple(3) do i
            padding_pixels(boundary_conditions[i], resolution[i])
        end

        # Calculate field attenuation pixels from boundary conditions
        # For AbsorbingBoundaryCondition, use attenuation_thickness
        # For PeriodicBoundaryCondition, use 0
        field_attenuation_pixel = ntuple(3) do i
            bc = boundary_conditions[i]
            if bc isa AbsorbingBoundaryCondition
                round(Int, bc.attenuation_thickness / (resolution[i] * 2))
            else
                0
            end
        end

        # Calculate region of interest
        ROI = (
            boundary_thickness_pixel[1] + 1, boundary_thickness_pixel[1] + grid_size[1],
            boundary_thickness_pixel[2] + 1, boundary_thickness_pixel[2] + grid_size[2],
            boundary_thickness_pixel[3] + 1, boundary_thickness_pixel[3] + grid_size[3]
        )

        new{T}(
            # Configuration fields
            wavelength, permittivity_bg, resolution, grid_size,
            boundary_conditions,  # NEW: Store boundary conditions directly
            iterations_max, tolerance, linear_solver, preconditioner,
            # Solver state fields
            backend,  # KernelAbstractions.jl backend (cached)
            nothing,  # Green_function - computed lazily (handles flip internally)
            nothing,  # potential - set during solve
            nothing,  # potential_device - set during solve
            nothing,  # permittivity - set during solve
            boundary_thickness_pixel,  # Derived from boundary conditions
            field_attenuation_pixel,   # Derived from boundary conditions
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
    ConvergentBornSolver(; kwargs...) -> ConvergentBornSolver

Public keyword-argument constructor for ConvergentBornSolver with boundary condition support.

# Required Arguments
- `wavelength::Real`: Wavelength in background medium (meters)
- `resolution::NTuple{3, <:Real}`: Spatial resolution (dx, dy, dz) in meters
- `grid_size::NTuple{3, Int}`: Grid dimensions (Nx, Ny, Nz)
- `boundary_conditions`: Either a single `AbstractBoundaryCondition` (applied to all dimensions)
  or `NTuple{3, AbstractBoundaryCondition}` for per-dimension control

# Optional Arguments
- `permittivity_bg::Real = 1.0`: Background relative permittivity
- `iterations_max::Int = -1`: Maximum Born iterations (-1 for auto)
- `tolerance::Real = 1e-6`: Convergence tolerance
- `linear_solver = KrylovJL_GMRES()`: LinearSolve.jl algorithm
- `preconditioner::Symbol = :none`: Preconditioning strategy (:none, :diagonal, :ilu)
- `device::Symbol = :cpu`: Compute device (:cpu, :cuda, :amdgpu, :metal, :oneapi)

# Examples
```julia
# Simple absorbing boundaries (same on all dimensions)
bc = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6,
    sharpness = 0.9
)
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32),
    boundary_conditions = bc
)

# Mixed boundaries: periodic in x,y, absorbing in z
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32),
    boundary_conditions = (
        PeriodicBoundaryCondition(),
        PeriodicBoundaryCondition(),
        bc
    )
)

# GPU acceleration with CUDA
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    boundary_conditions = bc,
    device = :cuda
)
```
"""
function ConvergentBornSolver(;
        wavelength::Real,
        permittivity_bg::Real = 1.0,
        resolution::NTuple{3, <:Real},
        grid_size::NTuple{3, Int},
        boundary_conditions,  # Accept any tuple or single boundary condition
        iterations_max::Int = -1,
        tolerance::Real = 1e-6,
        linear_solver::SciMLLinearSolveAlgorithm = KrylovJL_GMRES(),
        preconditioner::Symbol = :none,
        device::Symbol = :cpu)
    # Determine floating-point type from inputs (promote to highest precision)
    T = promote_type(
        typeof(float(wavelength)),
        typeof(float(permittivity_bg)),
        typeof(float(first(resolution))),
        typeof(float(tolerance))
    )

    # Convert boundary_conditions to NTuple{3, AbstractBoundaryCondition{T}}
    bcs = if boundary_conditions isa AbstractBoundaryCondition
        # Single boundary condition - apply to all dimensions
        # Convert to correct type parameter
        bc_converted = _convert_boundary_type(boundary_conditions, T)
        (bc_converted, bc_converted, bc_converted)
    else
        # Per-dimension boundary conditions - convert each
        ntuple(3) do i
            _convert_boundary_type(boundary_conditions[i], T)
        end
    end

    # Create backend from device symbol
    backend = create_backend_from_symbol(device)

    # Call internal constructor with converted types
    return ConvergentBornSolver{T}(
        T(wavelength),
        T(permittivity_bg),
        T.(resolution),
        grid_size,
        bcs,
        iterations_max,
        T(tolerance),
        linear_solver,
        preconditioner,
        backend
    )
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

The CBS iteration `Efield_{n+1} = Efield_0 + A * Efield_n` can be reformulated as
a linear system `(1-G*V*) * Efield = Efield_0` where:
- G*V* = (G + G_flip)/2 * V
- V is the potential (permittivity contrast)
- G, G_flip are dyadic Green's functions

This operator represents only the core CBS operator without the preconditioner R.
The preconditioner R = (1i/ε_imag) * potential is applied separately through LinearSolve.jl.

# Fields
- `solver::ConvergentBornSolver{T}`: Reference to the CBS solver
- `temp_array::AbstractArray`: Preallocated temporary 4D field array for efficiency
"""
struct CBSLinearOperator{T <: AbstractFloat}
    solver::ConvergentBornSolver{T}
    temp_array::AbstractArray  # Preallocated temporary 4D field array

    function CBSLinearOperator(solver::ConvergentBornSolver{T}) where {T}
        # Preallocate temporary array to match field dimensions
        field_size = size(solver.potential)
        temp_array = KernelAbstractions.zeros(solver.backend, Complex{T}, field_size..., 3)

        new{T}(solver, temp_array)
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
    psi = op.temp_array

    # Reshape flattened input to 4D field array
    field_shape = size(psi)
    x_field = reshape(x, field_shape)

    # Apply Green's operator: G*V*x = (G + G_flip)/2 * V * x

    # Step 1: V * x (element-wise multiplication with potential)
    psi .= potential .* x_field

    # Step 2: (G + G_flip)/2
    conv!(solver.Green_function, psi)

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
Base.size(op::CBSLinearOperator) = (length(op.temp_array), length(op.temp_array))
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
- `EMfield::ElectromagneticField`: Electromagnetic field structure containing E and H fields

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
EMfield = solve(solver, source, permittivity)

# Multi-source (coherent interference)
sources = [source1, source2, source3]
EMfield = solve(solver, sources, permittivity)
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

    # Initialize solver state (computes potential V with proper padding)
    _initialize_solver!(solver, permittivity)

    # Generate incident field
    incident_field = generate_incident_field(sources, solver)

    # Generate source term: V * E_incident + i*eps_imag * E_incident
    source_term = solver.potential .* incident_field.E .+
                  Complex{T}(0, solver.eps_imag) .* incident_field.E

    # Transfer source term to device
    source_term = to_device(solver.backend, source_term)

    # Solve scattering equation using CBS algorithm
    scattered_field = _solve_cbs_scattering(solver, source_term)

    # Transfer to host, add incident field, and crop to ROI
    total_field = to_host(scattered_field) + incident_field
    result_field = crop_to_ROI(total_field, solver)

    return result_field
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
    k_bg_medium = T(2π) * sqrt(solver.permittivity_bg) / solver.wavelength
    potential_unpadded = Complex{T}.(k_bg_medium^2 *
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

    # Transfer potential to device
    solver.potential_device = to_device(solver.backend, solver.potential)

    # Calculate optimal iteration count if automatic
    if solver.iterations_max < 0
        steps = abs(2 * k_bg_medium / solver.eps_imag)
        domain_extent = norm(size(solver.potential)[1:3] .* solver.resolution)
        Bornmax_opt = ceil(Int, domain_extent / steps / 2 + 1) * 2
        solver.Bornmax = Bornmax_opt * abs(solver.iterations_max)
    else
        solver.Bornmax = solver.iterations_max
    end

    # Initialize Green's function
    k_square = k_bg_medium^2 + Complex{T}(0, solver.eps_imag)
    solver.Green_function = DyadicGreen(
        solver.potential_device, k_square, solver.resolution, solver.boundary_conditions)

    # Reset iteration tracking
    solver.iteration_count = 0
    empty!(solver.residual_history)

    return nothing
end

"""
    _solve_cbs_scattering(solver, source) -> ElectromagneticField

Solve the electromagnetic scattering problem using LinearSolve.jl CBS implementation.

This function reformulates the Convergent Born Series as a linear system (I - A) * Efield = Efield_0
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
- `ElectromagneticField{T}`: Scattered electromagnetic field (device-resident)
"""
function _solve_cbs_scattering(
        solver::ConvergentBornSolver{T},
        source::AbstractArray{Complex{T}, 4}
) where {T <: AbstractFloat}

    # Prepare right-hand side: Field_0 from first CBS iteration
    size_field = size(source)
    Green_fn = solver.Green_function
    psi = KernelAbstractions.zeros(solver.backend, Complex{T}, size_field)

    if solver.Bornmax >= 1
        # Compute Field_0: first iteration of CBS without R preconditioner
        # Field_0 = (G + G_flip)/2 * source
        psi .= source

        # Apply Green's function
        conv!(Green_fn, psi)

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
        Efield = reshape(sol.u, size_field)

        # Update solver statistics if available
        if hasfield(typeof(sol), :iters)
            solver.iteration_count = sol.iters
        end
        if hasfield(typeof(sol), :resid)
            push!(solver.residual_history, sol.resid)
        end
    end

    # Compute magnetic field using same method as iterative version
    Hfield = _compute_magnetic_field(solver, Efield)

    # Get padded grid size from potential
    padded_grid_size = size(solver.potential)

    # Return ElectromagneticField (device-resident)
    return ElectromagneticField(Efield, Hfield, padded_grid_size, solver.resolution, solver.wavelength)
end

"""
    _compute_magnetic_field(solver, Efield) -> AbstractArray

Compute magnetic field from electric field using Maxwell's equation:
H = -i/k₀ * (n₀/Z₀) * ∇ × E

where k₀ is the free-space wave number and Z₀ = 377 Ω is the impedance.
"""
function _compute_magnetic_field(
        solver::ConvergentBornSolver{T},
        Efield::AbstractArray{Complex{T}, 4}
) where {T <: AbstractFloat}

    # Create curl operator for field with padding
    array_type = typeof(solver.potential)
    curl_op = Curl(array_type, T, size(solver.potential), solver.resolution)

    # Compute curl of E field
    Hfield = conv(curl_op, Efield)

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

Apply field attenuation masks to the solver's potential dimension-by-dimension.

Each boundary condition applies attenuation along its corresponding dimension.
This allows mixed boundary types (e.g., periodic in x,y and absorbing in z).
"""
function _apply_attenuation_to_potential!(solver::ConvergentBornSolver{T}) where {T}
    # Apply attenuation dimension-by-dimension using boundary condition logic
    for dim in 1:3
        bc = solver.boundary_conditions[dim]

        # Skip dimensions with no padding/attenuation
        if !requires_padding(bc)
            continue
        end

        # For absorbing boundaries, apply attenuation along this dimension
        if bc isa AbsorbingBoundaryCondition
            max_L = solver.boundary_thickness_pixel[dim]
            L = min(max_L, solver.field_attenuation_pixel[dim])

            if max_L == 0
                continue
            end

            # Create attenuation window using the boundary condition's profile
            window = _create_attenuation_window(bc.profile, L, T)

            # Apply sharpness factor from boundary condition
            window = window .* bc.sharpness .+ (1 - bc.sharpness)

            # Create full filter
            roi_size = solver.ROI[2 * dim] - solver.ROI[2 * dim - 1] + 1
            remaining_padding = max_L - L

            filter_1d = vcat(
                window,
                ones(T, roi_size + 2*remaining_padding),
                reverse(window)
            )

            # Apply mask directly to potential along the specified dimension
            if dim == 1
                solver.potential .*= reshape(filter_1d, :, 1, 1)
            elseif dim == 2
                solver.potential .*= reshape(filter_1d, 1, :, 1)
            else  # dim == 3
                solver.potential .*= reshape(filter_1d, 1, 1, :)
            end
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

    # Display boundary conditions
    bc_x = solver.boundary_conditions[1]
    bc_y = solver.boundary_conditions[2]
    bc_z = solver.boundary_conditions[3]
    if bc_x == bc_y == bc_z
        # Same boundary condition on all dimensions
        print(io, "\n  boundary_conditions: $(typeof(bc_x).name.name) (all dimensions)")
    else
        # Different boundary conditions per dimension
        print(io, "\n  boundary_conditions:")
        print(io, "\n    x: $(typeof(bc_x).name.name)")
        print(io, "\n    y: $(typeof(bc_y).name.name)")
        print(io, "\n    z: $(typeof(bc_z).name.name)")
    end

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
