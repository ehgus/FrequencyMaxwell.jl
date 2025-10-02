"""
Field generation utilities for electromagnetic sources.

This module provides utilities for generating incident electromagnetic fields
from sources, including multi-source superposition, and validation.
"""

"""
    generate_incident_field(sources::Vector{<:AbstractCurrentSource}, solver) -> ElectromagneticField

Generate incident electromagnetic field from multiple coherent sources with proper superposition.

This method implements coherent field superposition where multiple electromagnetic sources
interfere in the same simulation domain. Complex amplitudes are summed to produce
combined incident fields with constructive and destructive interference.

Each source generates its field on the full padded grid, ensuring physically correct
incident fields in the boundary regions.

# Arguments
- `solver`: Solver object containing grid configuration and boundary conditions
- `sources::Vector{<:AbstractCurrentSource{T}}`: Vector of coherent electromagnetic sources

# Returns
- `ElectromagneticField{T}`: Total electromagnetic field with interference on padded grid

# Algorithm
1. Validate source compatibility (same wavelength for coherent interference)
2. Generate incident field from first source (on padded grid)
3. Coherently add fields from remaining sources: E_total = Σ E_i, H_total = Σ H_i

This produces proper electromagnetic interference patterns where sources interfere
constructively (bright fringes) or destructively (dark fringes).

# Example
```julia
sources = [source1, source2, source3]
incident_field = generate_incident_field(sources, solver)
E_incident = incident_field.E  # Already on padded grid with interference
```
"""
function generate_incident_field(
        sources::Vector{<:AbstractCurrentSource{T}},
        solver
) where {T <: AbstractFloat}
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

    # Generate field from first source (establishes grid size)
    total_field = generate_incident_field(sources[1], solver)

    # Sum contributions from remaining sources (coherent superposition)
    for i in 2:length(sources)
        # Generate field for this source (already on padded grid)
        incident_field = generate_incident_field(sources[i], solver)

        # Coherent addition using ElectromagneticField + operator
        total_field = total_field + incident_field
    end

    return total_field
end
