"""
Plane wave electromagnetic source implementation.

This module provides plane wave sources for electromagnetic simulations,
including polarized plane waves with arbitrary propagation directions.
"""

"""
    PlaneWaveSource{T<:AbstractFloat} <: AbstractCurrentSource{T}

Electromagnetic plane wave source with specified polarization and propagation direction.

This struct provides a type-safe, validated plane wave source that generates
incident electromagnetic fields satisfying Maxwell's equations.

# Type Parameters
- `T<:AbstractFloat`: Floating-point precision type

# Fields
- `wavelength::T`: Wavelength in the background medium (meters)
- `polarization::SVector{3, Complex{T}}`: Complex polarization vector (V/m)
- `k_vector::SVector{3, T}`: Wave vector direction (normalized)
- `amplitude::T`: Field amplitude (V/m)
- `phase::T`: Initial phase offset (radians)

# Constructor
```julia
PlaneWaveSource(;
    wavelength::Real,
    polarization::AbstractVector{<:Number},
    k_vector::AbstractVector{<:Real},
    amplitude::Real = 1.0,
    phase::Real = 0.0
)
```

# Validation
The constructor automatically validates:
- `wavelength > 0`
- `polarization ⟂ k_vector` (transverse field condition)
- `norm(k_vector) ≈ 1` (normalized wave vector)
- `amplitude ≥ 0`

# Example
```julia
# X-polarized plane wave propagating in +Z direction
source = PlaneWaveSource(
    wavelength = 500e-9,           # 500 nm
    polarization = [1.0, 0.0, 0.0], # Linear X polarization
    k_vector = [0.0, 0.0, 1.0],     # +Z propagation
    amplitude = 1.0                  # 1 V/m amplitude
)

# Circularly polarized wave
source_circular = PlaneWaveSource(
    wavelength = 633e-9,
    polarization = [1.0, 1.0im, 0.0]/√2,  # Right circular polarization
    k_vector = [0.0, 0.0, 1.0],
    amplitude = 1.0
)
```
"""
struct PlaneWaveSource{T <: AbstractFloat} <: AbstractCurrentSource{T}
    wavelength::T
    polarization::SVector{3, Complex{T}}
    k_vector::SVector{3, T}
    amplitude::T
    phase::T

    function PlaneWaveSource{T}(
            wavelength::T,
            polarization::SVector{3, Complex{T}},
            k_vector::SVector{3, T},
            amplitude::T,
            phase::T
    ) where {T <: AbstractFloat}

        # Validate parameters
        wavelength > 0 || throw(ArgumentError("wavelength must be positive"))
        amplitude ≥ 0 || throw(ArgumentError("amplitude must be non-negative"))

        # Check wave vector normalization
        k_norm = norm(k_vector)
        abs(k_norm - 1) < 1e-10 ||
            throw(ArgumentError("k_vector must be normalized (|k| = 1)"))

        # Check transverse field condition (E ⟂ k)
        E_dot_k = real(dot(polarization, k_vector))
        abs(E_dot_k) < 1e-10 ||
            throw(ArgumentError("polarization must be transverse to k_vector"))

        return new{T}(wavelength, polarization, k_vector, amplitude, phase)
    end
end

# Convenience constructor with keyword arguments and type promotion
function PlaneWaveSource(;
        wavelength::Real,
        polarization::AbstractVector{<:Number},
        k_vector::AbstractVector{<:Real},
        amplitude::Real = 1.0,
        phase::Real = 0.0
)
    # Promote types
    T = promote_type(typeof(wavelength), typeof(amplitude), typeof(phase),
        real(eltype(polarization)), eltype(k_vector))
    T <: AbstractFloat || (T = Float64)

    # Convert to static vectors with proper types
    pol_static = SVector{3, Complex{T}}(Complex{T}.(polarization))
    k_static = SVector{3, T}(T.(k_vector))

    return PlaneWaveSource{T}(
        T(wavelength), pol_static, k_static, T(amplitude), T(phase)
    )
end

# Implement required interface methods

"""
    source_wavelength(source::PlaneWaveSource) -> T

Get the wavelength of the plane wave source.
"""
source_wavelength(source::PlaneWaveSource) = source.wavelength

"""
    cross(a::SVector{3}, b::SVector{3}) -> SVector{3}

Cross product for 3D static vectors.
"""
function cross(a::SVector{3, T1}, b::SVector{3, T2}) where {T1, T2}
    T = promote_type(T1, T2)
    return SVector{3, T}(
        a[2]*b[3] - a[3]*b[2],
        a[3]*b[1] - a[1]*b[3],
        a[1]*b[2] - a[2]*b[1]
    )
end

"""
    show(io::IO, source::PlaneWaveSource)

Custom display for plane wave source objects.
"""
function Base.show(io::IO, source::PlaneWaveSource{T}) where {T}
    print(io, "PlaneWaveSource{$T}:")
    print(io, "\n  wavelength: $(source.wavelength)")
    print(io, "\n  polarization: $(source.polarization)")
    print(io, "\n  k_vector: $(source.k_vector)")
    print(io, "\n  amplitude: $(source.amplitude)")
    print(io, "\n  phase: $(source.phase)")
end

"""
    generate_incident_field(source::PlaneWaveSource, solver) -> ElectromagneticField

Generate incident plane wave field on the solver's padded computational grid.

This function generates electromagnetic plane wave field satisfying Maxwell's equations:
- E(r) = E₀ * exp(ik·r + iφ)
- H(r) = (k × E) / (ωμ₀) = (k × E) / Z₀k₀

The field is generated on the full padded grid (including boundary regions) because
a plane wave physically exists everywhere in space, not just in the computational ROI.
The phase reference center remains at the center of the original (unpadded) domain
to ensure consistent field phases.

# Arguments
- `source::PlaneWaveSource{T}`: Plane wave source configuration
- `solver`: Solver object containing grid configuration (grid_size, resolution, permittivity_bg)

# Returns
- `ElectromagneticField{T}`: Electromagnetic field on padded grid

# Example
```julia
incident_field = generate_incident_field(plane_wave_source, solver)
E_array = incident_field.E  # Access electric field array (padded)
```
"""
function generate_incident_field(
        source::PlaneWaveSource{T},
        solver
) where {T <: AbstractFloat}

    # Extract grid parameters from solver
    grid_size = solver.grid_size  # Original computational grid
    resolution = solver.resolution
    permittivity_bg = solver.permittivity_bg
    padding = ntuple(3) do i
        padding_pixels(solver.boundary_conditions[i], resolution[i])
    end

    # Compute padded grid size
    padded_grid_size = grid_size .+ 2 .* padding

    # Initialize field arrays on padded grid
    E_incident = zeros(Complex{T}, padded_grid_size..., 3)
    H_incident = zeros(Complex{T}, padded_grid_size..., 3)

    # Wave parameters
    k0 = T(2π / source.wavelength)  # Free-space wave number
    k_bg = k0 * sqrt(permittivity_bg)  # Background wave number

    # Normalize k-vector and polarization
    k_hat = source.k_vector ./ norm(source.k_vector)  # Unit k-vector
    k_vec = k_bg .* k_hat  # Full k-vector in medium
    E0 = source.polarization ./ norm(source.polarization)  # Normalized polarization

    # Ensure E ⟂ k (transversality condition)
    E_perp = E0 .- (dot(E0, k_hat) * k_hat)  # Remove parallel component
    E_perp = E_perp ./ norm(E_perp)  # Renormalize

    # Compute H from E using H = k × E / (ωμ₀) = k × E / (k₀Z₀)
    H_perp = cross(k_hat, E_perp) ./ T(377 / sqrt(permittivity_bg))

    # Calculate domain center for phase reference (center of original unpadded domain)
    center = ntuple(i -> (grid_size[i] - 1) * resolution[i] / 2, 3)

    # Compute plane wave fields for all grid points (including padded regions)
    for i in 1:padded_grid_size[1], j in 1:padded_grid_size[2], k in 1:padded_grid_size[3]
        # Position relative to domain center (accounting for padding offset)
        r = [
            (i - 1 - padding[1]) * resolution[1] - center[1],
            (j - 1 - padding[2]) * resolution[2] - center[2],
            (k - 1 - padding[3]) * resolution[3] - center[3]
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

    # Create ElectromagneticField object
    return ElectromagneticField(
        E_incident,
        H_incident,
        padded_grid_size,
        resolution,
        source.wavelength
    )
end
