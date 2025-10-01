"""
Boundary condition abstraction for electromagnetic solvers.

This module provides a clean separation between boundary condition logic and the
core electromagnetic solver implementation, enabling flexible boundary condition
configuration and future extensions.

# Design Philosophy

Boundary conditions in electromagnetic simulations control:
1. **Padding requirements** - how much extra space is needed around the computational domain
2. **Green's function parameters** - subpixel shifts for proper boundary treatment
3. **Field attenuation** - masks that prevent boundary reflections

The Convergent Born Series (CBS) algorithm always uses Green's function averaging
`(G + G_flip)/2` regardless of boundary type. Boundary conditions control the
**subpixel shift values** (0.0 for periodic, ±0.25 for absorbing), not the
averaging itself.

# Interface

All boundary conditions must implement:
- `padding_pixels(bc, resolution)` - calculate padding needed in pixels
- `subpixel_shift(bc)` - return Green's function subpixel shift
- `requires_padding(bc)` - whether this boundary needs array padding
- `apply_attenuation!(bc, potential, params)` - apply field attenuation masks

# Example

```julia
# Create boundary conditions
periodic = PeriodicBoundaryCondition{Float64}()
absorbing = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,           # 2 micron padding
    attenuation_thickness = 1.5e-6,  # 1.5 micron attenuation
    sharpness = 0.9               # Sharp transition
)

# Query properties
shift_periodic = subpixel_shift(periodic)      # 0.0
shift_absorbing = subpixel_shift(absorbing)    # 0.25
padding = padding_pixels(absorbing, 50e-9)     # 40 pixels
```
"""

# Abstract Type Hierarchy

"""
    AbstractBoundaryCondition{T<:AbstractFloat}

Abstract supertype for all electromagnetic boundary conditions.

Boundary conditions in electromagnetic simulations control how fields behave at
the computational domain boundaries. Different boundary types provide different
physical behaviors:

- **Periodic**: Fields wrap around (useful for infinite periodic structures)
- **Absorbing**: Fields decay to zero (prevents reflections from boundaries)
- **PML**: Perfectly matched layers (future implementation)
- **PEC/PMC**: Perfect electric/magnetic conductors (future implementation)

# Type Parameter
- `T<:AbstractFloat`: Floating-point precision type (Float32, Float64, etc.)

# Interface Requirements

All concrete boundary condition types must implement:

## Required Methods

### Query Methods (read-only, no side effects):
- `padding_pixels(bc::AbstractBoundaryCondition, resolution::T) -> Int`
  Calculate the number of padding pixels needed for this boundary condition.

- `subpixel_shift(bc::AbstractBoundaryCondition) -> T`
  Return the subpixel shift for Green's function creation.
  This affects how the dyadic Green's function handles boundaries.

- `requires_padding(bc::AbstractBoundaryCondition) -> Bool`
  Whether this boundary condition requires array padding.

### Action Methods (may modify state):
- `apply_attenuation!(bc::AbstractBoundaryCondition, potential, params...) -> Nothing`
  Apply field attenuation masks to the potential array.
  For boundaries without attenuation (e.g., periodic), this is a no-op.

# Implementation Notes

The Convergent Born Series (CBS) algorithm uses Green's function averaging:
```
ψ_scattered = (G + G_flip)/2 * V * ψ
```
where G_flip uses opposite subpixel shift from G. This averaging is intrinsic
to CBS and happens **regardless** of boundary type. Boundary conditions control
the shift values (0.0 for periodic, ±0.25 for absorbing), not the averaging.

# See Also
- `PeriodicBoundaryCondition`: Periodic boundary implementation
- `AbsorbingBoundaryCondition`: Absorbing boundary implementation
"""
abstract type AbstractBoundaryCondition{T <: AbstractFloat} end

# Periodic Boundary Condition

"""
    PeriodicBoundaryCondition{T} <: AbstractBoundaryCondition{T}

Periodic boundary condition for electromagnetic simulations.

Periodic boundaries treat the computational domain as infinite and repeating.
Fields at one boundary edge are identical to fields at the opposite edge,
enabling simulation of infinite periodic structures like gratings, photonic
crystals, and metamaterials.

# Mathematical Description

For a periodic boundary in the x-direction:
```
E(x=0, y, z) = E(x=L, y, z)
E(x+L, y, z) = E(x, y, z)  for all x
```

# Green's Function Treatment

Periodic boundaries use **zero subpixel shift** in the dyadic Green's function:
```
subpixel_shift = 0.0
```

The CBS algorithm still uses Green's function averaging `(G + G_flip)/2`, but
both G and G_flip have the same shift (zero), making them identical for periodic
boundaries.

# Characteristics
- **Padding**: No additional padding required (but current implementation maintains compatibility)
- **Subpixel shift**: 0.0 (no shift for periodic boundaries)
- **Attenuation**: None (fields are not attenuated)
- **Physical meaning**: Simulates infinite periodic structure

# Constructor

```julia
PeriodicBoundaryCondition{Float64}()  # Explicit type
PeriodicBoundaryCondition{T}() where T  # Generic
```

# Example

```julia
# Create periodic boundary
bc = PeriodicBoundaryCondition{Float64}()

# Query properties
shift = subpixel_shift(bc)      # Returns 0.0
needs_pad = requires_padding(bc)  # Returns false
padding = padding_pixels(bc, 50e-9)  # Returns 0

# Apply attenuation (no-op for periodic)
apply_attenuation!(bc, potential, params...)  # Does nothing
```

# See Also
- `AbsorbingBoundaryCondition`: For non-periodic boundaries with field absorption
- `AbstractBoundaryCondition`: Base interface documentation
"""
struct PeriodicBoundaryCondition{T <: AbstractFloat} <: AbstractBoundaryCondition{T}
    # Minimal marker type - no fields needed
    # All behavior is defined through multiple dispatch on interface methods
end

"""
    PeriodicBoundaryCondition(::Type{T}=Float64) -> PeriodicBoundaryCondition{T}

Convenience constructor for PeriodicBoundaryCondition with optional type specification.

# Arguments
- `T::Type{<:AbstractFloat} = Float64`: Floating-point precision type (optional)

# Examples
```julia
# Default Float64 precision
bc = PeriodicBoundaryCondition()

# Explicit Float32 precision
bc32 = PeriodicBoundaryCondition(Float32)
```
"""
function PeriodicBoundaryCondition(::Type{T} = Float64) where {T <: AbstractFloat}
    PeriodicBoundaryCondition{T}()
end

# Attenuation Profile Types

"""
    AttenuationProfile

Enumeration of available attenuation profile functions for absorbing boundaries.

# Available Profiles

- **TanhProfile**: Hyperbolic tangent profile (default)
  - Smooth, continuous attenuation with excellent stability
  - Formula: `window(x) = [tanh(x) - tanh(-2.5)] / 2`
  - Best for: General-purpose electromagnetic simulations

- **ExponentialProfile**: Exponential decay profile
  - Rapid decay with natural physical interpretation
  - Formula: `window(x) = exp(-((L-x)/w)²)`
  - Best for: Strong absorption, wide-band simulations

- **TukeyProfile**: Tukey (tapered cosine) window
  - Cosine-based smooth transition
  - Formula: `window(x) = 0.5 * (1 + cos(π*(1-x/L)))`
  - Best for: Minimal dispersion, narrowband problems

# Example

```julia
# Create absorbing boundary with different profiles
bc_tanh = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6,
    profile = TanhProfile  # Default
)

bc_exp = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6,
    profile = ExponentialProfile  # Rapid absorption
)
```

# References
- Tukey window: J. W. Tukey, "An introduction to the calculations of numerical spectrum analysis"
- Planck-taper: McKechan et al., "A tapering window for time-domain templates"
"""
@enum AttenuationProfile begin
    TanhProfile          # Hyperbolic tangent (default, current implementation)
    ExponentialProfile   # Exponential decay
    TukeyProfile         # Tukey (tapered cosine) window
end

# Absorbing Boundary Condition

"""
    AbsorbingBoundaryCondition{T} <: AbstractBoundaryCondition{T}

Absorbing boundary condition with configurable attenuation profiles.

Absorbing boundaries prevent reflections from computational domain edges by
smoothly attenuating electromagnetic fields to zero near boundaries. Multiple
attenuation profiles are available (Tanh, Exponential, Tukey, Planck-taper),
each optimized for different simulation scenarios.

# Mathematical Description

Field attenuation is applied dimension-by-dimension using the selected profile.
The general form is:
```
attenuation(x) = sharpness * window(x) + (1 - sharpness)
```

where:
- x ∈ [0, L] is position within attenuation layer
- window(x) is the profile-specific attenuation function
- sharpness ∈ [0, 1] controls transition abruptness

See `AttenuationProfile` documentation for specific profile formulas.

# Green's Function Treatment

Absorbing boundaries use **non-zero subpixel shift** for proper boundary treatment:
```
subpixel_shift = 0.25
flip_subpixel_shift = -0.25
```

The CBS algorithm averages `(G + G_flip)/2` where these shifts ensure proper
handling of non-periodic boundary conditions.

# Fields
- `thickness::T`: Physical thickness of boundary padding in meters
- `attenuation_thickness::T`: Physical thickness of attenuation layer in meters
- `sharpness::T`: Sharpness factor ∈ [0, 1] controlling attenuation profile
  - 0.0: Gentle, gradual attenuation (minimal artifacts, wider transition)
  - 1.0: Sharp, abrupt attenuation (narrow transition, possible artifacts)
  - 0.7-0.9: Recommended range for most applications
- `profile::AttenuationProfile`: Attenuation function type (default: TanhProfile)

# Constructor

```julia
AbsorbingBoundaryCondition(;
    thickness::Real,
    attenuation_thickness::Real,
    sharpness::Real = 0.9,
    profile::AttenuationProfile = TanhProfile
)
```

# Validation

Constructor validates:
- All parameters must be positive: `thickness > 0`, `attenuation_thickness > 0`
- Attenuation must fit within padding: `attenuation_thickness ≤ thickness`
- Sharpness must be in valid range: `0 ≤ sharpness ≤ 1`

# Example

```julia
# Create absorbing boundary with default Tanh profile
bc = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,           # 2 micron total padding
    attenuation_thickness = 1.5e-6,  # 1.5 micron attenuation layer
    sharpness = 0.9               # Sharp attenuation (90%)
)

# Use different attenuation profile
bc_exp = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6,
    sharpness = 0.85,
    profile = ExponentialProfile  # Exponential decay
)

# Query properties
shift = subpixel_shift(bc)           # Returns 0.25
needs_pad = requires_padding(bc)     # Returns true
padding = padding_pixels(bc, 50e-9)  # Returns 40 pixels (2μm / 50nm)

# Apply attenuation to potential
apply_attenuation!(bc, potential, boundary_thickness_pixel,
                  field_attenuation_pixel, ROI)
```

# Performance Considerations

- Attenuation is applied dimension-by-dimension for efficiency
- Masks are computed on-the-fly and immediately applied (no storage overhead)
- GPU-compatible through broadcasting operations
- Tanh computation is vectorized for performance

# See Also
- `PeriodicBoundaryCondition`: For periodic boundaries without attenuation
- `AbstractBoundaryCondition`: Base interface documentation
- `apply_attenuation!`: Detailed attenuation algorithm documentation
"""
struct AbsorbingBoundaryCondition{T <: AbstractFloat} <: AbstractBoundaryCondition{T}
    thickness::T                    # Physical thickness in meters
    attenuation_thickness::T        # Attenuation layer thickness in meters
    sharpness::T                    # Sharpness factor (0-1)
    profile::AttenuationProfile     # Attenuation function type

    function AbsorbingBoundaryCondition{T}(
            thickness::T,
            attenuation_thickness::T,
            sharpness::T,
            profile::AttenuationProfile
    ) where {T <: AbstractFloat}
        # Validate parameters
        thickness > 0 ||
            throw(ArgumentError("thickness must be positive, got: $thickness"))
        attenuation_thickness > 0 ||
            throw(ArgumentError("attenuation_thickness must be positive, got: $attenuation_thickness"))
        attenuation_thickness ≤ thickness ||
            throw(ArgumentError(
                "attenuation_thickness ($attenuation_thickness) must be ≤ thickness ($thickness)"
            ))
        0 ≤ sharpness ≤ 1 ||
            throw(ArgumentError("sharpness must be in [0, 1], got: $sharpness"))

        new{T}(thickness, attenuation_thickness, sharpness, profile)
    end
end

"""
    AbsorbingBoundaryCondition(; kwargs...) -> AbsorbingBoundaryCondition

Convenience constructor with keyword arguments and automatic type promotion.

# Arguments
- `thickness::Real`: Physical thickness of boundary padding (meters)
- `attenuation_thickness::Real`: Physical thickness of attenuation layer (meters)
- `sharpness::Real = 0.9`: Sharpness factor ∈ [0, 1]
- `profile::AttenuationProfile = TanhProfile`: Attenuation function type

# Example
```julia
bc = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6,
    sharpness = 0.9,
    profile = TanhProfile
)
```
"""
function AbsorbingBoundaryCondition(;
        thickness::Real,
        attenuation_thickness::Real,
        sharpness::Real = 0.9,
        profile::AttenuationProfile = TanhProfile
)
    # Promote to common floating-point type
    T = promote_type(typeof(thickness), typeof(attenuation_thickness), typeof(sharpness))
    T <: AbstractFloat || (T = Float64)

    return AbsorbingBoundaryCondition{T}(T(thickness), T(attenuation_thickness), T(sharpness), profile)
end

# Interface Implementation - Query Methods

"""
    padding_pixels(bc::AbstractBoundaryCondition, resolution::Real) -> Int

Calculate the number of padding pixels needed for a boundary condition.

Padding pixels are additional grid points added around the computational domain
to properly implement boundary conditions. The padding size depends on the
boundary type and the spatial resolution.

# Arguments
- `bc::AbstractBoundaryCondition`: Boundary condition instance
- `resolution::Real`: Spatial resolution (grid spacing) in meters

# Returns
- `Int`: Number of padding pixels required in one direction

# Method Implementations

## PeriodicBoundaryCondition
Returns 0 pixels - periodic boundaries don't require padding (fields wrap around).
Note: Current implementation maintains compatibility with existing code that may
still pad periodic boundaries.

## AbsorbingBoundaryCondition
Returns pixels needed for the specified thickness:
```julia
padding_pixels = round(Int, thickness / (resolution * 2))
```
The factor of 2 accounts for both sides of the boundary layer.

# Example
```julia
resolution = 50e-9  # 50 nm grid spacing

# Periodic boundary
bc_periodic = PeriodicBoundaryCondition{Float64}()
pixels_p = padding_pixels(bc_periodic, resolution)  # Returns 0

# Absorbing boundary with 2 μm thickness
bc_absorb = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6
)
pixels_a = padding_pixels(bc_absorb, resolution)  # Returns 20 pixels
# Calculation: 2.0e-6 / (50e-9 * 2) = 20
```

# See Also
- `requires_padding`: Check if padding is needed
- `subpixel_shift`: Get Green's function shift parameter
"""
function padding_pixels(bc::PeriodicBoundaryCondition{T}, resolution::Real) where {T}
    # Periodic boundaries don't require padding (fields wrap around)
    return 0
end

function padding_pixels(bc::AbsorbingBoundaryCondition{T}, resolution::Real) where {T}
    # Calculate padding pixels from physical thickness
    # Factor of 2 accounts for both sides of boundary
    return round(Int, bc.thickness / (T(resolution) * 2))
end

"""
    subpixel_shift(bc::AbstractBoundaryCondition) -> T

Return the subpixel shift for dyadic Green's function creation.

The subpixel shift parameter controls how the dyadic Green's function handles
boundaries in the electromagnetic scattering calculation. Different boundary
types require different shifts for proper physical behavior.

# Physical Meaning

The subpixel shift affects phase ramps applied during Green's function convolution:
```
phase_ramp = exp(i * k · shift)
```

This shift ensures proper boundary condition enforcement in the CBS algorithm's
Green's function averaging: `(G + G_flip)/2` where G_flip uses opposite shift.

# Returns
- `T<:AbstractFloat`: Subpixel shift value in units of grid spacing

# Method Implementations

## PeriodicBoundaryCondition
Returns `0.0` - no shift needed for periodic boundaries.
Both G and G_flip have zero shift, making them identical.

## AbsorbingBoundaryCondition
Returns `0.25` - quarter-pixel shift for non-periodic boundaries.
The flipped Green's function uses `-0.25`, and averaging reduces artifacts.

# Example
```julia
# Periodic boundary
bc_p = PeriodicBoundaryCondition{Float64}()
shift_p = subpixel_shift(bc_p)  # Returns 0.0

# Absorbing boundary
bc_a = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6
)
shift_a = subpixel_shift(bc_a)  # Returns 0.25

# In CBS solver, this becomes:
# Green_fn = DyadicGreen(..., shifts = shift_a)      # +0.25
# flip_Green_fn = DyadicGreen(..., shifts = -shift_a) # -0.25
# result = (conv(Green_fn) + conv(flip_Green_fn)) / 2
```

# Implementation Note

The CBS algorithm ALWAYS uses Green's function averaging `(G + G_flip)/2`
regardless of boundary type. This method only controls the shift values, not
whether averaging occurs.

# See Also
- `DyadicGreen`: Dyadic Green's function implementation using subpixel shifts
- `padding_pixels`: Calculate required padding for boundary condition
"""
function subpixel_shift(bc::PeriodicBoundaryCondition{T}) where {T}
    # Periodic boundaries use zero shift
    return T(0)
end

function subpixel_shift(bc::AbsorbingBoundaryCondition{T}) where {T}
    # Non-periodic boundaries use quarter-pixel shift
    # CBS averages G(+0.25) and G(-0.25) for proper boundary treatment
    return T(0.25)
end

"""
    requires_padding(bc::AbstractBoundaryCondition) -> Bool

Check whether a boundary condition requires array padding.

Some boundary conditions need extra grid points around the computational domain
(padding), while others don't. This method indicates whether the array size
should be extended to accommodate the boundary condition.

# Returns
- `Bool`: `true` if padding is required, `false` otherwise

# Method Implementations

## PeriodicBoundaryCondition
Returns `false` - periodic boundaries wrap around and don't need padding.
However, current implementation may still pad for compatibility.

## AbsorbingBoundaryCondition
Returns `true` - absorbing boundaries need padding for attenuation layers.

# Example
```julia
# Periodic boundary
bc_p = PeriodicBoundaryCondition{Float64}()
requires_padding(bc_p)  # Returns false

# Absorbing boundary
bc_a = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6
)
requires_padding(bc_a)  # Returns true
```

# See Also
- `padding_pixels`: Calculate the actual padding size
"""
function requires_padding(bc::PeriodicBoundaryCondition)
    return false
end

function requires_padding(bc::AbsorbingBoundaryCondition)
    return true
end

"""
    requires_averaging(bc::AbstractBoundaryCondition) -> Bool

Check whether a boundary condition requires Green's function averaging.

The Convergent Born Series uses Green's function averaging `(G + G_flip)/2`
where G and G_flip use opposite subpixel shifts. This method indicates whether
the boundary condition requires this averaging mechanism.

# Returns
- `Bool`: `true` if averaging is required, `false` otherwise

# Method Implementations

## PeriodicBoundaryCondition
Returns `false` - periodic boundaries use zero shift, making G = G_flip.
Averaging still occurs in implementation but has no effect.

## AbsorbingBoundaryCondition
Returns `true` - non-periodic boundaries require averaging of ±0.25 shifts.

# Example
```julia
bc_p = PeriodicBoundaryCondition{Float64}()
requires_averaging(bc_p)  # Returns false

bc_a = AbsorbingBoundaryCondition(thickness = 2.0e-6, attenuation_thickness = 1.5e-6)
requires_averaging(bc_a)  # Returns true
```

# Implementation Note
This is used by `DyadicGreen` to determine whether to internally compute
the flip Green's function and perform averaging.
"""
function requires_averaging(bc::PeriodicBoundaryCondition)
    return false
end

function requires_averaging(bc::AbsorbingBoundaryCondition)
    return true
end

# Profile-Specific Attenuation Window Functions

"""
    _create_attenuation_window(profile::AttenuationProfile, L::Int, T::Type) -> Vector{T}

Create an attenuation window for the specified profile type.

# Arguments
- `profile::AttenuationProfile`: Attenuation function type
- `L::Int`: Number of points in attenuation layer
- `T::Type{<:AbstractFloat}`: Floating-point type

# Returns
- `Vector{T}`: Attenuation window values from 0 (fully attenuated) to 1 (no attenuation)
"""
function _create_attenuation_window(profile::AttenuationProfile, L::Int, ::Type{T}) where {T <:
                                                                                           AbstractFloat}
    if L <= 0
        return T[]
    end

    if profile == TanhProfile
        # Hyperbolic tangent profile (original implementation)
        tanh_vals = tanh.(range(T(-2.5), T(2.5), length = L))
        return (tanh_vals ./ tanh(T(2.5)) .- tanh(T(-2.5))) ./ 2

    elseif profile == ExponentialProfile
        # Exponential decay: exp(-((L-x)/w)²)
        x = range(T(0), T(1), length = L)
        # Use Gaussian decay with width parameter
        w = T(0.3)  # Width parameter controls decay rate
        return exp.(-(((T(1) .- x) ./ w) .^ 2))

    elseif profile == TukeyProfile
        # Tukey window (tapered cosine): 0.5 * (1 + cos(π*(1-x/L)))
        x = range(T(0), T(1), length = L)
        return T(0.5) .* (1 .+ cos.(π .* (1 .- x)))

    else
        error("Unknown attenuation profile: $profile")
    end
end

# Interface Implementation - Action Methods

"""
    apply_attenuation!(
        bc::AbstractBoundaryCondition,
        potential::AbstractArray{<:Number, 3},
        boundary_thickness_pixel::NTuple{3, Int},
        field_attenuation_pixel::NTuple{3, Int},
        ROI::NTuple{6, Int}
    ) -> Nothing

Apply field attenuation masks to the potential array for boundary condition enforcement.

Field attenuation prevents electromagnetic field reflections from computational
domain boundaries by smoothly reducing field amplitudes near edges. This method
modifies the potential array in-place, applying dimension-specific attenuation
profiles.

# Arguments
- `bc::AbstractBoundaryCondition`: Boundary condition instance
- `potential::AbstractArray{<:Number, 3}`: Potential array V (will be modified in-place)
- `boundary_thickness_pixel::NTuple{3, Int}`: Boundary padding in pixels per dimension
- `field_attenuation_pixel::NTuple{3, Int}`: Attenuation layer thickness in pixels per dimension
- `ROI::NTuple{6, Int}`: Region of interest bounds (x_min, x_max, y_min, y_max, z_min, z_max)

# Method Implementations

## PeriodicBoundaryCondition
No-op (does nothing). Periodic boundaries don't attenuate fields.

## AbsorbingBoundaryCondition
Applies tanh-based attenuation masks dimension-by-dimension:

### Algorithm (per dimension):
1. Calculate attenuation layer size: `L = min(max_L, attenuation_thickness_pixel)`
2. Create tanh window: `window = [tanh(x) - tanh(-2.5)] / 2` for x ∈ [-2.5, 2.5]
3. Apply sharpness: `window = sharpness * window + (1 - sharpness)`
4. Create full 1D filter: `[window, ones(ROI_size), reverse(window)]`
5. Broadcast multiply along dimension: `potential .*= filter_1d`

### Attenuation Profile:
```
          ┌─────────────────────────┐
          │   Region of Interest    │
    ┌─────┼─────────────────────────┼─────┐
    │Atten│                         │Atten│
    └─────┴─────────────────────────┴─────┘
    ↑                                      ↑
    0.0                                    1.0
    (fully attenuated)                     (no attenuation)
```

# Example
```julia
# Setup
potential = Complex{Float64}.(randn(100, 100, 100))
boundary_thickness_pixel = (20, 20, 20)
field_attenuation_pixel = (15, 15, 15)
ROI = (21, 80, 21, 80, 21, 80)

# Periodic boundary - no change
bc_p = PeriodicBoundaryCondition{Float64}()
apply_attenuation!(bc_p, potential, boundary_thickness_pixel,
                  field_attenuation_pixel, ROI)
# potential is unchanged

# Absorbing boundary - applies attenuation
bc_a = AbsorbingBoundaryCondition(
    thickness = 2.0e-6,
    attenuation_thickness = 1.5e-6,
    sharpness = 0.9
)
apply_attenuation!(bc_a, potential, boundary_thickness_pixel,
                  field_attenuation_pixel, ROI)
# potential now has smooth attenuation near boundaries
```

# Performance Notes
- Operates dimension-by-dimension for memory efficiency
- No intermediate mask arrays stored (computed and applied directly)
- GPU-compatible through broadcasting
- Tanh computation is vectorized

# See Also
- `AbsorbingBoundaryCondition`: Boundary type with attenuation
- `PeriodicBoundaryCondition`: Boundary type without attenuation
"""
function apply_attenuation!(
        bc::PeriodicBoundaryCondition,
        potential::AbstractArray{<:Number, 3},
        boundary_thickness_pixel::NTuple{3, Int},
        field_attenuation_pixel::NTuple{3, Int},
        ROI::NTuple{6, Int}
)
    # Periodic boundaries don't attenuate fields - this is a no-op
    return nothing
end

function apply_attenuation!(
        bc::AbsorbingBoundaryCondition{T},
        potential::AbstractArray{<:Number, 3},
        boundary_thickness_pixel::NTuple{3, Int},
        field_attenuation_pixel::NTuple{3, Int},
        ROI::NTuple{6, Int}
) where {T}
    # Apply attenuation dimension-by-dimension
    for dim in 1:3
        max_L = boundary_thickness_pixel[dim]
        L = min(max_L, field_attenuation_pixel[dim])

        if max_L == 0
            continue
        end

        # Create attenuation window using selected profile
        window = _create_attenuation_window(bc.profile, L, T)

        # Apply sharpness factor
        # sharpness = 1.0: use full window (sharp transition)
        # sharpness = 0.0: use all ones (no attenuation)
        window = window .* bc.sharpness .+ (1 - bc.sharpness)

        # Create full 1D filter with attenuation on both ends
        roi_size = ROI[2 * dim] - ROI[2 * dim - 1] + 1
        remaining_padding = max_L - L

        filter_1d = vcat(
            window,                              # Left attenuation
            ones(T, roi_size + 2 * remaining_padding),  # ROI (no attenuation)
            reverse(window)                      # Right attenuation
        )

        # Apply 1D filter along the specified dimension via broadcasting
        if dim == 1
            potential .*= reshape(filter_1d, :, 1, 1)
        elseif dim == 2
            potential .*= reshape(filter_1d, 1, :, 1)
        else  # dim == 3
            potential .*= reshape(filter_1d, 1, 1, :)
        end
    end

    return nothing
end

# Type Conversion Utilities

"""
    _convert_boundary_type(bc::AbstractBoundaryCondition, ::Type{T}) -> AbstractBoundaryCondition{T}

Internal: Convert boundary condition to the target floating-point type.

This utility function ensures type consistency when creating solvers with
mixed-precision inputs. It converts boundary condition parameters to the
target floating-point type while preserving all other properties.

# Arguments
- `bc::AbstractBoundaryCondition`: Boundary condition to convert
- `T::Type{<:AbstractFloat}`: Target floating-point type

# Returns
- Boundary condition with parameters converted to type `T`

# Examples
```julia
# Convert Float32 boundary condition to Float64
bc32 = AbsorbingBoundaryCondition{Float32}(2.0f-6, 1.5f-6, 0.9f0, TanhProfile)
bc64 = _convert_boundary_type(bc32, Float64)
```
"""
function _convert_boundary_type(bc::PeriodicBoundaryCondition, ::Type{T}) where {T <:
                                                                                 AbstractFloat}
    return PeriodicBoundaryCondition{T}()
end

function _convert_boundary_type(bc::AbsorbingBoundaryCondition, ::Type{T}) where {T <:
                                                                                  AbstractFloat}
    return AbsorbingBoundaryCondition{T}(
        T(bc.thickness),
        T(bc.attenuation_thickness),
        T(bc.sharpness),
        bc.profile
    )
end

# Display Methods (REPL-friendly)

"""
    show(io::IO, bc::PeriodicBoundaryCondition)

Display periodic boundary condition in REPL-friendly format.
"""
function Base.show(io::IO, bc::PeriodicBoundaryCondition{T}) where {T}
    print(io, "PeriodicBoundaryCondition{$T}()")
end

"""
    show(io::IO, bc::AbsorbingBoundaryCondition)

Display absorbing boundary condition with parameter details.
"""
function Base.show(io::IO, bc::AbsorbingBoundaryCondition{T}) where {T}
    print(io, "AbsorbingBoundaryCondition{$T}:")
    print(io, "\n  thickness: $(bc.thickness) m")
    print(io, "\n  attenuation_thickness: $(bc.attenuation_thickness) m")
    print(io, "\n  sharpness: $(bc.sharpness)")
    print(io, "\n  profile: $(bc.profile)")
    print(io, "\n  subpixel_shift: $(subpixel_shift(bc))")
end
