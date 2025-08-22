"""
Abstract source types and interfaces for electromagnetic field generation.

This module defines the abstract interface that all electromagnetic sources
must implement, providing a consistent API for field generation.
"""

"""
    generate_incident_fields(source::AbstractCurrentSource, grid_config) -> (E_field, H_field)

Generate incident electromagnetic fields for the given source on a computational grid.

# Arguments
- `source::AbstractCurrentSource{T}`: Source object containing source parameters
- `grid_config`: Grid configuration object containing spatial discretization

# Returns
- `E_field::AbstractArray{Complex{T}, 4}`: Electric field (last dimension = 3 for components)
- `H_field::AbstractArray{Complex{T}, 4}`: Magnetic field (last dimension = 3 for components)

# Requirements
All concrete source types must implement this method to generate physically
consistent electromagnetic fields satisfying Maxwell's equations.
"""
function generate_incident_fields end

"""
    source_power(source::AbstractCurrentSource) -> Real

Calculate the total power delivered by the electromagnetic source.

This is an optional interface that sources can implement for power normalization
and energy conservation checks.
"""
function source_power end

"""
    source_wavelength(source::AbstractCurrentSource) -> Real

Get the primary wavelength of the electromagnetic source.

All sources must have a well-defined wavelength for the electromagnetic solver.
"""
function source_wavelength end

"""
    validate_source(source::AbstractCurrentSource) -> Bool

Validate that the source parameters are physically consistent.

This function checks for common issues such as:
- Non-physical parameter values
- Inconsistent polarization and propagation directions
- Invalid wavelength or frequency values

Returns `true` if validation passes, throws an `ArgumentError` if validation fails.
"""
function validate_source end
