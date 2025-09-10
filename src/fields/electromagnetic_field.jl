"""
Electromagnetic field data structures and utilities.

This module provides structured representations for electromagnetic fields
with associated grid information and physical units.
"""

"""
    ElectromagneticField{T<:AbstractFloat, N, AT<:AbstractArray{Complex{T}, N}}

Structured container for electromagnetic field data with grid information.

This struct provides a type-safe wrapper for electromagnetic field arrays,
ensuring proper dimensionality and providing convenient access to field
components and derived quantities.

# Type Parameters
- `T<:AbstractFloat`: Floating-point precision type
- `N`: Array dimensionality (typically 4 for 3D fields + 3 components)
- `AT<:AbstractArray`: Array type (for CPU/GPU flexibility)

# Fields
- `E::AT`: Electric field array (V/m)
- `H::AT`: Magnetic field array (A/m)
- `grid_size::NTuple{3, Int}`: Spatial grid dimensions
- `resolution::NTuple{3, T}`: Spatial resolution (meters)
- `wavelength::T`: Wavelength in the medium (meters)

# Constructor
```julia
ElectromagneticField(E_array, H_array, grid_size, resolution, wavelength)
```

# Example
```julia
# Create field data structure
E_field = zeros(Complex{Float64}, 128, 128, 32, 3)  # 3D grid + 3 components
H_field = zeros(Complex{Float64}, 128, 128, 32, 3)
grid_size = (128, 128, 32)
resolution = (50e-9, 50e-9, 100e-9)  # 50nm x 50nm x 100nm voxels
wavelength = 500e-9  # 500 nm

fields = ElectromagneticField(E_field, H_field, grid_size, resolution, wavelength)

# Access field components
Ex = fields.E[:, :, :, 1]  # X-component of electric field
Hy = fields.H[:, :, :, 2]  # Y-component of magnetic field
```
"""
struct ElectromagneticField{T <: AbstractFloat, N, AT <: AbstractArray{Complex{T}, N}}
    E::AT
    H::AT
    grid_size::NTuple{3, Int}
    resolution::NTuple{3, T}
    wavelength::T

    function ElectromagneticField{T, N, AT}(
            E::AT,
            H::AT,
            grid_size::NTuple{3, Int},
            resolution::NTuple{3, T},
            wavelength::T
    ) where {T <: AbstractFloat, N, AT <: AbstractArray{Complex{T}, N}}

        # Validate field array dimensions
        size(E) == size(H) || throw(ArgumentError("E and H arrays must have same size"))

        # For 4D arrays (3D space + 3 components), check grid consistency
        if N == 4
            size(E)[1:3] == grid_size ||
                throw(ArgumentError("field array spatial dimensions must match grid_size"))
            size(E)[4] == 3 ||
                throw(ArgumentError("field arrays must have 3 components in last dimension"))
        end

        # Validate physical parameters
        all(grid_size .> 0) ||
            throw(ArgumentError("grid_size must have positive dimensions"))
        all(resolution .> 0) || throw(ArgumentError("resolution must be positive"))
        wavelength > 0 || throw(ArgumentError("wavelength must be positive"))

        return new{T, N, AT}(E, H, grid_size, resolution, wavelength)
    end
end

# Convenience constructor with automatic type inference
function ElectromagneticField(
        E::AbstractArray{Complex{T}, N},
        H::AbstractArray{Complex{T}, N},
        grid_size::NTuple{3, Int},
        resolution::NTuple{3, <:Real},
        wavelength::Real
) where {T <: AbstractFloat, N}

    # Promote resolution and wavelength to match field precision
    resolution_T = NTuple{3, T}(T.(resolution))
    wavelength_T = T(wavelength)

    AT = typeof(E)
    return ElectromagneticField{T, N, AT}(E, H, grid_size, resolution_T, wavelength_T)
end

"""
    domain_size(fields::ElectromagneticField) -> NTuple{3, T}

Calculate the total physical domain size.
"""
function domain_size(fields::ElectromagneticField{T}) where {T}
    return ntuple(i -> fields.grid_size[i] * fields.resolution[i], 3)
end

"""
    field_energy(fields::ElectromagneticField) -> T

Calculate the total electromagnetic energy density in the domain.

The energy density is given by:
u = (1/2) * (ε₀|E|² + μ₀|H|²)

For simplicity, this assumes vacuum permittivity and permeability.
"""
function field_energy(fields::ElectromagneticField{T}) where {T}
    ε₀ = T(8.854187817e-12)  # Vacuum permittivity (F/m)
    μ₀ = T(4π * 1e-7)        # Vacuum permeability (H/m)

    # Calculate |E|² and |H|² summed over components
    E_energy = sum(abs2, fields.E)
    H_energy = sum(abs2, fields.H)

    # Total energy density
    energy_density = (ε₀ * E_energy + μ₀ * H_energy) / 2

    # Multiply by voxel volume
    voxel_volume = prod(fields.resolution)

    return energy_density * voxel_volume
end

"""
    poynting_vector(fields::ElectromagneticField) -> AbstractArray{Complex{T}, 4}

Calculate the complex Poynting vector S = E × H*.

Returns an array with the same spatial dimensions as the input fields,
with the last dimension representing the 3 vector components.
"""
function poynting_vector(fields::ElectromagneticField{T, 4}) where {T}
    # Get field components
    Ex, Ey, Ez = fields.E[:, :, :, 1], fields.E[:, :, :, 2], fields.E[:, :, :, 3]
    Hx, Hy, Hz = fields.H[:, :, :, 1], fields.H[:, :, :, 2], fields.H[:, :, :, 3]

    # Calculate S = E × H* (complex conjugate of H)
    Sx = Ey .* conj.(Hz) .- Ez .* conj.(Hy)
    Sy = Ez .* conj.(Hx) .- Ex .* conj.(Hz)
    Sz = Ex .* conj.(Hy) .- Ey .* conj.(Hx)

    # Combine into 4D array
    S = similar(fields.E)
    S[:, :, :, 1] = Sx
    S[:, :, :, 2] = Sy
    S[:, :, :, 3] = Sz

    return S
end

"""
    field_intensity(fields::ElectromagneticField; component=:total) -> AbstractArray{T, 3}

Calculate the field intensity |E|² or individual component intensities.

# Arguments
- `fields::ElectromagneticField`: Input electromagnetic fields
- `component::Symbol`: Component to calculate (:total, :x, :y, :z)

# Returns
- `intensity::AbstractArray{T, 3}`: Field intensity distribution
"""
function field_intensity(fields::ElectromagneticField{T, 4}; component::Symbol = :total) where {T}
    if component == :total
        return sum(abs2, fields.E; dims = 4)[:, :, :, 1]
    elseif component == :x
        return abs2.(fields.E[:, :, :, 1])
    elseif component == :y
        return abs2.(fields.E[:, :, :, 2])
    elseif component == :z
        return abs2.(fields.E[:, :, :, 3])
    else
        throw(ArgumentError("component must be :total, :x, :y, or :z"))
    end
end

"""
    extract_plane(fields::ElectromagneticField, plane_axis::Int, plane_index::Int) -> ElectromagneticField

Extract a 2D plane from the 3D electromagnetic field.

# Arguments
- `fields::ElectromagneticField`: Input 3D fields
- `plane_axis::Int`: Axis normal to the plane (1=X, 2=Y, 3=Z)
- `plane_index::Int`: Index along the plane axis

# Returns
- `plane_fields::ElectromagneticField`: 2D electromagnetic field on the specified plane
"""
function extract_plane(
        fields::ElectromagneticField{T, 4},
        plane_axis::Int,
        plane_index::Int
) where {T}
    1 ≤ plane_axis ≤ 3 || throw(ArgumentError("plane_axis must be 1, 2, or 3"))
    1 ≤ plane_index ≤ fields.grid_size[plane_axis] ||
        throw(ArgumentError("plane_index out of bounds"))

    # Extract plane data
    if plane_axis == 1  # YZ plane
        E_plane = fields.E[plane_index:plane_index, :, :, :]
        H_plane = fields.H[plane_index:plane_index, :, :, :]
        plane_grid_size = (1, fields.grid_size[2], fields.grid_size[3])
        plane_resolution = (
            fields.resolution[1], fields.resolution[2], fields.resolution[3])
    elseif plane_axis == 2  # XZ plane
        E_plane = fields.E[:, plane_index:plane_index, :, :]
        H_plane = fields.H[:, plane_index:plane_index, :, :]
        plane_grid_size = (fields.grid_size[1], 1, fields.grid_size[3])
        plane_resolution = (
            fields.resolution[1], fields.resolution[2], fields.resolution[3])
    else  # XY plane
        E_plane = fields.E[:, :, plane_index:plane_index, :]
        H_plane = fields.H[:, :, plane_index:plane_index, :]
        plane_grid_size = (fields.grid_size[1], fields.grid_size[2], 1)
        plane_resolution = (
            fields.resolution[1], fields.resolution[2], fields.resolution[3])
    end

    return ElectromagneticField(
        E_plane, H_plane, plane_grid_size, plane_resolution, fields.wavelength)
end

"""
    show(io::IO, fields::ElectromagneticField)

Custom display for electromagnetic field objects.
"""
function Base.show(io::IO, fields::ElectromagneticField{T, N, AT}) where {T, N, AT}
    print(io, "ElectromagneticField{$T, $N, $(AT.name)}:")
    print(io, "\n  grid_size: $(fields.grid_size)")
    print(io, "\n  resolution: $(fields.resolution)")
    print(io, "\n  domain_size: $(domain_size(fields))")
    print(io, "\n  wavelength: $(fields.wavelength)")
    print(io, "\n  field_energy: $(field_energy(fields))")
end
