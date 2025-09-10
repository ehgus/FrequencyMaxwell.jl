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

# Configure solver
config = ConvergentBornConfig(
    wavelength = 500e-9,
    permittivity_bg = 1.33^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 32)
)

# Create solver
solver = ConvergentBornSolver(config)

# Define source
source = PlaneWaveSource(
    wavelength = config.wavelength,
    polarization = [1.0, 0.0, 0.0],
    k_vector = [0.0, 0.0, 1.0]
)

# Solve electromagnetic problem
E_field, H_field = solve(solver, source, permittivity_distribution)
```
"""
module FrequencyMaxwell

using LinearAlgebra
using LinearSolve
using FFTW
using StaticArrays
using KernelAbstractions
using SciMLBase

# Export core types and functions
export
# Configuration types
      ConvergentBornConfig,

# Solver types
      ConvergentBornSolver,

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
      generate_incident_fields

# Include submodules in correct order (dependencies first)
include("core/types.jl")
include("core/configuration.jl")
include("core/dyadic_green.jl")
include("core/curl.jl")
include("sources/abstract_source.jl")
include("sources/plane_wave.jl")
include("fields/electromagnetic_field.jl")
include("geometry/phantoms.jl")
include("solvers/convergent_born.jl")

end # module FrequencyMaxwell
