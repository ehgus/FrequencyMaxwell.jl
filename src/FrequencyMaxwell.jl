"""
    FrequencyMaxwell

Modern Julia package for electromagnetic field simulation using frequency-domain Maxwell equations.
This package provides automatic differentiation-enabled solvers for electromagnetic scattering,
inverse design, and optimization problems.

# Key Features
- Type-safe, immutable configuration with mutable solver state
- LinearSolve.jl integration for flexible linear algebra backends
- Native automatic differentiation support via Zygote.jl
- Memory-efficient implementation with GPU acceleration support
- Comprehensive electromagnetic solver ecosystem

# Main Modules
- `Solvers`: Core electromagnetic solvers (ConvergentBorn, etc.)
- `Sources`: Electromagnetic source definitions (plane waves, focused beams)
- `Materials`: Material property database and modeling
- `Geometry`: Phantom generation and geometric utilities
- `Optimization`: Inverse design and optimization frameworks

# Example Usage
```julia
using FrequencyMaxwell

# Create solver with streamlined API (recommended)
solver = ConvergentBornSolver(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32)
)

# Define source
source = PlaneWaveSource(
    wavelength = solver.config.wavelength,
    polarization = [1.0, 0.0, 0.0],
    k_vector = [0.0, 0.0, 1.0]
)

# Solve electromagnetic problem
Efield, Hfield = solve(solver, source, permittivity_distribution)
```
"""
module FrequencyMaxwell

using LinearAlgebra
using LinearSolve
using FFTW
using StaticArrays
using KernelAbstractions
using SciMLBase
using Reexport
using Requires

# Reexport LinearSolve symbols to avoid name clashes
@reexport using LinearSolve: solve

# Export core types and functions
export
# Solver types
      ConvergentBornSolver,
      reset!,

# Source types  
      PlaneWaveSource,
      AbstractCurrentSource,

# Core functions
      solve,

# Phantom generation
      phantom_bead,
      phantom_plate,
      phantom_cylinder,

# Field utilities
      ElectromagneticField,
      field_energy,
      field_intensity,
      poynting_vector,
      extract_plane,

# Configuration utilities
      domain_size,
      grid_spacing,
      wavenumber_background,

# Source utilities
      source_wavelength,
      source_power,
      validate_source,
      generate_incident_fields,

# GPU Backend utilities
      is_backend_available

# Include submodules in correct order (dependencies first)
include("core/types.jl")
include("core/dyadic_green.jl")
include("core/curl.jl")
include("sources/abstract_source.jl")
include("sources/plane_wave.jl")
include("fields/electromagnetic_field.jl")
include("geometry/phantoms.jl")
include("solvers/convergent_born.jl")

# GPU Backend Registry for conditional loading
const GPU_BACKENDS = Dict{Symbol, Function}()

"""
    register_gpu_backend!(device::Symbol, constructor::Function)

Register a GPU backend constructor function for conditional loading.

This function allows GPU packages to register their backend constructors
when they are loaded via Requires.jl @require macros.
"""
function register_gpu_backend!(device::Symbol, constructor::Function)
    GPU_BACKENDS[device] = constructor
    return nothing
end

"""
    is_backend_available(device::Symbol) -> Bool

Check if a GPU backend is available for the specified device.

Returns true if the corresponding GPU package has been loaded and
registered its backend constructor.
"""
function is_backend_available(device::Symbol)
    return haskey(GPU_BACKENDS, device)
end

# Conditional GPU backend loading
function __init__()
    @require CUDA="052768ef-5323-5732-b1bb-66c8b64840ba" begin
        register_gpu_backend!(:cuda, () -> CUDA.CUDABackend())
    end

    @require AMDGPU="21141c5a-9bdb-4563-92ae-f87d6854732e" begin
        register_gpu_backend!(:amdgpu, () -> AMDGPU.ROCBackend())
    end

    @require Metal="dde4c033-4e86-420c-a63e-0dd931031962" begin
        register_gpu_backend!(:metal, () -> Metal.MetalBackend())
    end

    @require oneAPI="8f75cd03-7ff8-4ecb-9b8f-daf728133b1b" begin
        register_gpu_backend!(:oneapi, () -> oneAPI.oneAPIBackend())
    end
end

end # module FrequencyMaxwell
