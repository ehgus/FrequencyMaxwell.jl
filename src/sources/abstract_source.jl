"""
Abstract source types and interfaces for electromagnetic field generation.

This module defines the abstract interface that all electromagnetic sources
must implement, providing a consistent API for field generation.
"""

"""
    generate_incident_field(source::AbstractCurrentSource, solver) -> ElectromagneticField

Generate incident electromagnetic field for the given source on the solver's computational grid.

This is the primary interface that all source types must implement. The method generates
incident electromagnetic field on the solver's unpadded grid. Padding for boundary
conditions is handled separately by field generation utilities.

# Arguments
- `source::AbstractCurrentSource{T}`: Source object containing source parameters
- `solver`: Solver object containing grid configuration (grid_size, resolution, permittivity_bg)

# Returns
- `ElectromagneticField{T}`: Electromagnetic field structure containing:
  - `E::AbstractArray{Complex{T}, 4}`: Electric field (last dimension = 3 for components)
  - `H::AbstractArray{Complex{T}, 4}`: Magnetic field (last dimension = 3 for components)
  - Grid information (size, resolution, wavelength)

# Requirements
All concrete source types must implement this method to generate physically
consistent electromagnetic fields satisfying Maxwell's equations:
1. Fields must be generated on the unpadded grid (solver.grid_size)
2. E and H must be properly related through Maxwell's equations
3. Fields must satisfy appropriate physical constraints (e.g., E ⟂ k for plane waves)
4. Return ElectromagneticField object with properly populated fields and grid info

# Example Implementation
```julia
function generate_incident_field(source::MySource{T}, solver) where {T}
    # Extract grid parameters
    grid_size = solver.grid_size
    resolution = solver.resolution

    # Generate field arrays
    E_field = zeros(Complex{T}, grid_size..., 3)
    H_field = zeros(Complex{T}, grid_size..., 3)

    # Fill arrays with source-specific physics
    # ... implementation ...

    # Return ElectromagneticField object
    return ElectromagneticField(E_field, H_field, resolution, source.wavelength)
end
```
"""
function generate_incident_field end

"""
    source_wavelength(source::AbstractCurrentSource) -> Real

Get the primary wavelength of the electromagnetic source.

All sources must have a well-defined wavelength for the electromagnetic solver.
"""
function source_wavelength end
