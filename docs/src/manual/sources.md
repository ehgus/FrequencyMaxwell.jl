# Electromagnetic Sources

FrequencyMaxwell.jl provides a flexible framework for defining electromagnetic sources. The package includes built-in source types and an extensible interface for creating custom sources.

## Source Interface

All electromagnetic sources implement the `AbstractCurrentSource` interface, which provides a consistent API for source definition and field generation.

### Core Interface Methods

Every electromagnetic source must implement:

```julia
# Get source wavelength
wavelength = source_wavelength(source)

# Get source power (for normalization)
power = source_power(source)

# Validate source parameters
is_valid = validate_source(source)

# Generate incident fields on a grid
E_field, H_field = generate_incident_fields(source, config)
```

This interface enables:
- **Consistent API**: All sources work with the same solver interface
- **Multiple dispatch**: Different source types handled automatically
- **Extensibility**: Easy to add new source types
- **Validation**: Automatic parameter checking

## PlaneWaveSource

The `PlaneWaveSource` represents a plane wave electromagnetic field with specified polarization and propagation direction.

### Basic Usage

```julia
source = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [1.0, 0.0, 0.0],  # x-polarized
    k_vector = [0.0, 0.0, 1.0],      # propagating in +z
    amplitude = 1.0,
    phase = 0.0
)
```

### Parameters

**`wavelength`** (required)
- **Type**: `T <: AbstractFloat`
- **Description**: Vacuum wavelength of the electromagnetic radiation
- **Units**: Meters
- **Example**: `500e-9` for 500 nm light

**`polarization`** (required)
- **Type**: `AbstractVector` (3 components)
- **Description**: Electric field polarization vector
- **Constraints**: Must be perpendicular to k_vector (`E ⊥ k`)
- **Normalization**: Automatically normalized to unit magnitude
- **Example**: `[1.0, 0.0, 0.0]` for x-polarization

**`k_vector`** (required)
- **Type**: `AbstractVector` (3 components)
- **Description**: Propagation direction vector
- **Constraints**: Automatically normalized to unit magnitude
- **Example**: `[0.0, 0.0, 1.0]` for +z propagation

**`amplitude`** (optional, default: `1.0`)
- **Type**: `T <: AbstractFloat`
- **Description**: Electric field amplitude
- **Units**: V/m (or dimensionless for normalized fields)
- **Example**: `1.0` for unit amplitude

**`phase`** (optional, default: `0.0`)
- **Type**: `T <: AbstractFloat`
- **Description**: Phase offset of the electromagnetic wave
- **Units**: Radians
- **Example**: `π/2` for 90° phase shift

### Polarization Types

#### Linear Polarization

```julia
# x-polarized
source_x = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [1.0, 0.0, 0.0],
    k_vector = [0.0, 0.0, 1.0]
)

# y-polarized
source_y = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [0.0, 1.0, 0.0],
    k_vector = [0.0, 0.0, 1.0]
)

# 45° linear polarization
source_45 = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [1.0, 1.0, 0.0],  # Automatically normalized
    k_vector = [0.0, 0.0, 1.0]
)
```

#### Circular Polarization

```julia
# Right circular polarization (when viewed along propagation)
source_rcp = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [1.0, 1im, 0.0],  # Complex polarization
    k_vector = [0.0, 0.0, 1.0]
)

# Left circular polarization
source_lcp = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [1.0, -1im, 0.0],
    k_vector = [0.0, 0.0, 1.0]
)
```

#### Elliptical Polarization

```julia
# Elliptical polarization with 2:1 axis ratio
source_elliptical = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [2.0, 1im, 0.0],
    k_vector = [0.0, 0.0, 1.0]
)
```

### Propagation Directions

#### Standard Directions

```julia
# +z propagation (most common)
k_z = [0.0, 0.0, 1.0]

# -z propagation (counter-propagating)
k_z_neg = [0.0, 0.0, -1.0]

# +x propagation
k_x = [1.0, 0.0, 0.0]

# +y propagation
k_y = [0.0, 1.0, 0.0]
```

#### Oblique Incidence

```julia
# 30° incidence in xz-plane
k_oblique = [sin(30° * π/180), 0.0, cos(30° * π/180)]

source_oblique = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [0.0, 1.0, 0.0],  # y-polarized (TE)
    k_vector = k_oblique
)

# For TM polarization, polarization must be in xz-plane
polarization_tm = [cos(30° * π/180), 0.0, -sin(30° * π/180)]
source_tm = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = polarization_tm,
    k_vector = k_oblique
)
```

### Field Generation

The `generate_incident_fields` function creates the incident electromagnetic fields:

```julia
config = ConvergentBornConfig(...)
source = PlaneWaveSource(...)

# Generate fields on the computational grid
E_inc, H_inc = generate_incident_fields(source, config)

# E_inc and H_inc are ElectromagneticField objects
println(size(E_inc.E))  # (3, nx, ny, nz) - vector field components
println(size(H_inc.H))  # (3, nx, ny, nz) - vector field components
```

The generated fields satisfy:
- **Maxwell's equations**: ∇ × E = -iωμH, ∇ × H = iωεE
- **Orthogonality**: E ⊥ H ⊥ k for plane waves
- **Proper normalization**: Based on source amplitude

### Validation and Error Handling

The `PlaneWaveSource` constructor automatically validates parameters:

```julia
# This will throw an error - E not perpendicular to k
try
    source = PlaneWaveSource(
        wavelength = 500e-9,
        polarization = [1.0, 0.0, 1.0],  # Has z-component
        k_vector = [0.0, 0.0, 1.0]       # z-propagation
    )
catch e
    println("Error: ", e)  # "Polarization must be perpendicular to k_vector"
end

# Manual validation
source = PlaneWaveSource(...)
if validate_source(source)
    println("Source is valid")
else
    println("Source validation failed")
end
```

## Multi-Source Simulations

FrequencyMaxwell.jl supports coherent superposition of multiple sources:

### Interference Patterns

```julia
# Two coherent sources with phase difference
source1 = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [1.0, 0.0, 0.0],
    k_vector = [0.0, 0.0, 1.0],
    phase = 0.0
)

source2 = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [1.0, 0.0, 0.0],
    k_vector = [0.0, 0.0, 1.0],
    phase = π  # 180° phase difference
)

# Solve with both sources
E_field, H_field = solve(solver, [source1, source2], phantom)
```

### Polarization Control

```julia
# Create circular polarization from linear components
source_x = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [1.0, 0.0, 0.0],
    k_vector = [0.0, 0.0, 1.0],
    phase = 0.0
)

source_y = PlaneWaveSource(
    wavelength = 500e-9,
    polarization = [0.0, 1.0, 0.0],
    k_vector = [0.0, 0.0, 1.0],
    phase = π/2  # 90° phase shift
)

# Results in right circular polarization
E_field, H_field = solve(solver, [source_x, source_y], phantom)
```

### Beam Shaping

```julia
# Create a focused beam by interfering plane waves
sources = PlaneWaveSource[]

# Add multiple plane waves with different angles
for angle in -10:2:10  # ±10° in 2° steps
    k_vec = [sin(angle * π/180), 0.0, cos(angle * π/180)]
    source = PlaneWaveSource(
        wavelength = 500e-9,
        polarization = [0.0, 1.0, 0.0],
        k_vector = k_vec,
        amplitude = cos(angle * π/180)^2  # Amplitude weighting
    )
    push!(sources, source)
end

# Solve with all sources
E_field, H_field = solve(solver, sources, phantom)
```

## Custom Source Types

You can create custom electromagnetic sources by implementing the `AbstractCurrentSource` interface:

### Interface Implementation

```julia
# Define your custom source type
struct MyCustomSource{T<:AbstractFloat} <: AbstractCurrentSource{T}
    wavelength::T
    # ... other parameters
end

# Implement required interface methods
FrequencyMaxwell.source_wavelength(source::MyCustomSource) = source.wavelength

FrequencyMaxwell.source_power(source::MyCustomSource) = 1.0  # Or compute actual power

function FrequencyMaxwell.validate_source(source::MyCustomSource)
    # Implement validation logic
    return source.wavelength > 0
end

function FrequencyMaxwell.generate_incident_fields(source::MyCustomSource, config::ConvergentBornConfig)
    # Generate E and H fields for your source
    # Return ElectromagneticField objects
    E_field = ElectromagneticField(E_array)
    H_field = ElectromagneticField(H_array)
    return E_field, H_field
end
```

### Example: Gaussian Beam Source

```julia
struct GaussianBeamSource{T<:AbstractFloat} <: AbstractCurrentSource{T}
    wavelength::T
    beam_waist::T
    polarization::Vector{Complex{T}}
    focus_position::Vector{T}
    k_vector::Vector{T}
end

function generate_incident_fields(source::GaussianBeamSource, config::ConvergentBornConfig)
    # Implement Gaussian beam field generation
    # This would involve computing the complex beam profile
    # and applying the appropriate phase and amplitude
    # ...
end
```

## Source Utilities

### Wavelength Conversion

```julia
# Convert between wavelength and frequency
λ = 500e-9  # meters
c = 3e8     # m/s
f = c / λ   # frequency in Hz
ω = 2π * f  # angular frequency in rad/s
```

### Refractive Index and Wavelength

```julia
# Wavelength in medium
λ_vacuum = 500e-9
n_medium = 1.33  # water
λ_medium = λ_vacuum / n_medium

# Always use vacuum wavelength for source definition
source = PlaneWaveSource(
    wavelength = λ_vacuum,  # Use vacuum wavelength
    # ... other parameters
)
```

### Power and Intensity

```julia
source = PlaneWaveSource(
    wavelength = 500e-9,
    amplitude = 1.0,  # V/m
    # ... other parameters
)

# Power calculation (in vacuum)
ε₀ = 8.854e-12  # F/m
c = 3e8         # m/s
power_density = 0.5 * ε₀ * c * abs(source.amplitude)^2  # W/m²
```

## Best Practices

### Source Definition

1. **Use vacuum wavelength**: Always specify wavelength in vacuum
2. **Validate orthogonality**: Ensure E ⊥ k for plane waves
3. **Consistent units**: Use SI units throughout (meters, seconds, etc.)
4. **Normalization**: Choose appropriate amplitude scaling

### Multi-Source Coherence

1. **Same wavelength**: Coherent sources must have identical wavelengths
2. **Phase relationships**: Carefully control relative phases
3. **Power normalization**: Consider total power when using multiple sources

### Performance Considerations

1. **Source complexity**: Complex sources require more computation
2. **Grid resolution**: Ensure adequate sampling of source features
3. **Boundary effects**: Account for source positioning relative to boundaries