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
- `polarization::ComplexVec3{T}`: Complex polarization vector (V/m)
- `k_vector::RealVec3{T}`: Wave vector direction (normalized)
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
    polarization::ComplexVec3{T}
    k_vector::RealVec3{T}
    amplitude::T
    phase::T

    function PlaneWaveSource{T}(
            wavelength::T,
            polarization::ComplexVec3{T},
            k_vector::RealVec3{T},
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
    pol_static = ComplexVec3{T}(Complex{T}.(polarization))
    k_static = RealVec3{T}(T.(k_vector))

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
    source_power(source::PlaneWaveSource) -> T

Calculate the power density of the plane wave (intensity).

For a plane wave, this returns the time-averaged Poynting vector magnitude:
S = (1/2) * |E|² / Z₀ where Z₀ ≈ 377 Ω is the impedance of free space.
"""
function source_power(source::PlaneWaveSource{T}) where {T}
    E_magnitude_squared = norm(source.polarization)^2 * source.amplitude^2
    Z0 = T(376.730313668)  # Impedance of free space (Ω)
    return E_magnitude_squared / (2 * Z0)
end

"""
    validate_source(source::PlaneWaveSource) -> Bool

Validate the plane wave source parameters.
"""
function validate_source(source::PlaneWaveSource{T}) where {T}
    # Check basic parameter validity
    source.wavelength > 0 || return false
    source.amplitude ≥ 0 || return false

    # Check wave vector normalization
    abs(norm(source.k_vector) - 1) < 1e-10 || return false

    # Check transverse condition
    E_dot_k = real(dot(source.polarization, source.k_vector))
    abs(E_dot_k) < 1e-10 || return false

    return true
end

"""
    generate_incident_fields(source::PlaneWaveSource, grid_config) -> (E_field, H_field)

Generate incident plane wave fields on a computational grid.

# Algorithm
For a plane wave with polarization E₀ and wave vector k:
- E(r) = E₀ * exp(ik·r + iφ)
- H(r) = (k × E) / (μ₀ω) = (k × E) / Z₀k₀

where Z₀ is the impedance of free space and k₀ = 2π/λ.
"""
function _generate_incident_fields_arrays(
        source::PlaneWaveSource{T},
        grid_config
) where {T <: AbstractFloat}

    # Extract grid parameters
    grid_size = grid_config.grid_size
    resolution = grid_config.resolution

    # Calculate wave number in background medium
    k0 = T(2π) * sqrt(grid_config.permittivity_bg) / source.wavelength
    Z0 = T(376.730313668)  # Impedance of free space

    # Initialize field arrays
    E_field = zeros(Complex{T}, grid_size..., 3)
    H_field = zeros(Complex{T}, grid_size..., 3)

    # Calculate phase array exp(ik·r)
    phase_array = _calculate_phase_array(source, grid_config, k0)

    # Calculate magnetic field polarization: H₀ = k × E₀ / Z₀
    H_polarization = cross(source.k_vector, source.polarization) / Z0

    # Fill field arrays
    for component in 1:3
        E_field[:, :, :, component] = source.amplitude * source.polarization[component] *
                                      phase_array
        H_field[:, :, :, component] = source.amplitude * H_polarization[component] *
                                      phase_array
    end

    return E_field, H_field
end

"""
    _calculate_phase_array(source::PlaneWaveSource, grid_config, k0) -> AbstractArray

Calculate the spatial phase array exp(ik·r + iφ) for the plane wave.
"""
function _calculate_phase_array(
        source::PlaneWaveSource{T},
        grid_config,
        k0::T
) where {T <: AbstractFloat}
    grid_size = grid_config.grid_size
    resolution = grid_config.resolution

    # Create coordinate arrays
    phase_array = zeros(Complex{T}, grid_size)

    # Calculate domain center for phase reference
    center = ntuple(i -> (grid_size[i] - 1) * resolution[i] / 2, 3)

    # Fill phase array
    for k in 1:grid_size[3], j in 1:grid_size[2], i in 1:grid_size[1]
        # Position relative to center
        r = SVector{3, T}(
            (i - 1) * resolution[1] - center[1],
            (j - 1) * resolution[2] - center[2],
            (k - 1) * resolution[3] - center[3]
        )

        # Phase: k₀(k̂·r) + φ
        phase = k0 * dot(source.k_vector, r) + source.phase
        phase_array[i, j, k] = exp(1im * phase)
    end

    return phase_array
end

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
    print(io, "\n  power density: $(source_power(source))")
end

"""
    generate_incident_fields(source::PlaneWaveSource, grid_config) -> ElectromagneticField

Generate incident plane wave fields and return as an ElectromagneticField object.

This is a convenience wrapper around the lower-level field generation function
that returns a properly structured ElectromagneticField object.
"""
function generate_incident_fields(
        source::PlaneWaveSource{T},
        grid_config
) where {T <: AbstractFloat}

    # Generate field arrays
    E_field, H_field = _generate_incident_fields_arrays(source, grid_config)

    # Create ElectromagneticField object
    return ElectromagneticField(
        E_field,
        H_field,
        grid_config.grid_size,
        grid_config.resolution,
        source.wavelength
    )
end
