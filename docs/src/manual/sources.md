# Electromagnetic Sources

FrequencyMaxwell.jl provides a flexible framework for defining electromagnetic sources. The package includes built-in source types and an extensible interface for creating custom sources.

## Source Interface

All electromagnetic sources implement the `AbstractCurrentSource` interface, which provides a consistent API for source definition and field generation.

### Core Interface Methods

Every electromagnetic source must implement:

```julia
# Get source wavelength
wavelength = source_wavelength(source)

# Generate incident field on the solver's grid (returns ElectromagneticField)
EMfield = generate_incident_field(source, solver)
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

The `generate_incident_field` function creates the incident electromagnetic fields:

```julia
solver = ConvergentBornSolver(...)
source = PlaneWaveSource(...)

# Generate fields on the solver's grid (including padding)
EMfield = generate_incident_field(source, solver)

# EMfield is an ElectromagneticField object
println(size(EMfield.E))  # (nx_padded, ny_padded, nz_padded, 3) - vector field components
println(size(EMfield.H))  # (nx_padded, ny_padded, nz_padded, 3) - vector field components
```

The generated fields satisfy:
- **Maxwell's equations**: ∇ × E = -iωμH, ∇ × H = iωεE
- **Orthogonality**: E ⊥ H ⊥ k for plane waves
- **Proper normalization**: Based on source amplitude
- **Physical correctness**: Fields generated on full padded grid (not zero-padded)

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
Efield, Hfield = solve(solver, [source1, source2], phantom)
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
Efield, Hfield = solve(solver, [source_x, source_y], phantom)
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
Efield, Hfield = solve(solver, sources, phantom)
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

function FrequencyMaxwell.generate_incident_field(source::MyCustomSource, solver)
    # Extract grid parameters from solver
    grid_size = solver.grid_size
    resolution = solver.resolution
    padding = ntuple(3) do i
        padding_pixels(solver.boundary_conditions[i], resolution[i])
    end
    padded_grid_size = grid_size .+ 2 .* padding

    # Generate E and H field arrays on padded grid
    # (Implementation depends on your specific source type)
    E_array = zeros(Complex{eltype(resolution)}, padded_grid_size..., 3)
    H_array = zeros(Complex{eltype(resolution)}, padded_grid_size..., 3)
    # ... fill arrays with your source-specific field computation ...

    # Return ElectromagneticField object
    return ElectromagneticField(E_array, H_array, padded_grid_size, resolution, source.wavelength)
end
```

### Example: Gaussian Beam Source

```julia
struct GaussianBeamSource{T<:AbstractFloat} <: AbstractCurrentSource{T}
    wavelength::T
    beam_waist::T
    polarization::SVector{3, Complex{T}}
    focus_position::SVector{3, T}
    k_vector::SVector{3, T}
end

FrequencyMaxwell.source_wavelength(source::GaussianBeamSource) = source.wavelength

function FrequencyMaxwell.generate_incident_field(source::GaussianBeamSource, solver)
    # Extract grid parameters
    grid_size = solver.grid_size
    resolution = solver.resolution
    padding = solver.boundary_thickness_pixel
    padded_grid_size = grid_size .+ 2 .* padding

    # Implement Gaussian beam field generation on padded grid
    # This would involve computing the complex beam profile
    # and applying the appropriate phase and amplitude
    # ...

    return ElectromagneticField(E_field, H_field, padded_grid_size, resolution, source.wavelength)
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