"""
Electromagnetic field data structures and utilities.

This module provides structured representations for electromagnetic fields
with associated grid information and physical units.
"""

"""
    ElectromagneticField{T<:AbstractFloat}

Structured container for electromagnetic field data with grid information.

This struct provides a type-safe wrapper for electromagnetic field arrays,
ensuring proper dimensionality and providing convenient access to field
components and derived quantities.

# Type Parameters
- `T<:AbstractFloat`: Floating-point precision type

# Fields
- `E::AbstractArray{Complex{T}}`: Electric field array (V/m)
- `H::AbstractArray{Complex{T}}`: Magnetic field array (A/m)
- `grid_size::NTuple{3, Int}`: Spatial grid dimensions
- `resolution::NTuple{3, T}`: Spatial resolution (meters)
- `wavelength::T`: Wavelength in the medium (meters)

# Constructor
```julia
ElectromagneticField(Efield, Hfield, grid_size, resolution, wavelength)
```

# Example
```julia
# Create field data structure
Efield = zeros(Complex{Float64}, 128, 128, 32, 3)  # 3D grid + 3 components
Hfield = zeros(Complex{Float64}, 128, 128, 32, 3)
grid_size = (128, 128, 32)
resolution = (50e-9, 50e-9, 100e-9)  # 50nm x 50nm x 100nm voxels
wavelength = 500e-9  # 500 nm

EMfield = ElectromagneticField(Efield, Hfield, grid_size, resolution, wavelength)

# Access field components
Ex = EMfield.E[:, :, :, 1]  # X-component of electric field
Hy = EMfield.H[:, :, :, 2]  # Y-component of magnetic field
```
"""
struct ElectromagneticField{T <: AbstractFloat}
    E::AbstractArray{Complex{T}}
    H::AbstractArray{Complex{T}}
    grid_size::NTuple{3, Int}
    resolution::NTuple{3, T}
    wavelength::T

    function ElectromagneticField{T}(
            Efield::AbstractArray{Complex{T}},
            Hfield::AbstractArray{Complex{T}},
            grid_size::NTuple{3, Int},
            resolution::NTuple{3, T},
            wavelength::T
    ) where {T <: AbstractFloat}

        # Validate field array dimensions
        size(Efield) == size(Hfield) ||
            throw(ArgumentError("Efield and Hfield arrays must have same size"))

        # For 4D arrays (3D space + 3 components), check grid consistency
        if ndims(Efield) == 4
            size(Efield)[1:3] == grid_size ||
                throw(ArgumentError("field array spatial dimensions must match grid_size"))
            size(Efield)[4] == 3 ||
                throw(ArgumentError("field arrays must have 3 components in last dimension"))
        end

        # Validate physical parameters
        all(grid_size .> 0) ||
            throw(ArgumentError("grid_size must have positive dimensions"))
        all(resolution .> 0) || throw(ArgumentError("resolution must be positive"))
        wavelength > 0 || throw(ArgumentError("wavelength must be positive"))

        return new{T}(Efield, Hfield, grid_size, resolution, wavelength)
    end
end

# Convenience constructor with automatic type inference
function ElectromagneticField(
        Efield::AbstractArray{Complex{T}},
        Hfield::AbstractArray{Complex{T}},
        grid_size::NTuple{3, Int},
        resolution::NTuple{3, <:Real},
        wavelength::Real
) where {T <: AbstractFloat}

    # Promote resolution and wavelength to match field precision
    resolution_T = NTuple{3, T}(T.(resolution))
    wavelength_T = T(wavelength)

    return ElectromagneticField{T}(Efield, Hfield, grid_size, resolution_T, wavelength_T)
end

"""
    domain_size(EMfield::ElectromagneticField) -> NTuple{3, T}

Calculate the total physical domain size.
"""
function domain_size(EMfield::ElectromagneticField{T}) where {T}
    return ntuple(i -> EMfield.grid_size[i] * EMfield.resolution[i], 3)
end

"""
    field_energy(EMfield::ElectromagneticField) -> T

Calculate the total electromagnetic energy density in the domain.

The energy density is given by:
u = (1/2) * (ε₀|E|² + μ₀|H|²)

For simplicity, this assumes vacuum permittivity and permeability.
"""
function field_energy(EMfield::ElectromagneticField{T}) where {T}
    ε₀ = T(8.854187817e-12)  # Vacuum permittivity (F/m)
    μ₀ = T(4π * 1e-7)        # Vacuum permeability (H/m)

    # Calculate |E|² and |H|² summed over components
    E_energy = sum(abs2, EMfield.E)
    H_energy = sum(abs2, EMfield.H)

    # Total energy density
    energy_density = (ε₀ * E_energy + μ₀ * H_energy) / 2

    # Multiply by voxel volume
    voxel_volume = prod(EMfield.resolution)

    return energy_density * voxel_volume
end

"""
    poynting_vector(EMfield::ElectromagneticField) -> AbstractArray{Complex{T}}

Calculate the complex Poynting vector S = E × H*.

Returns a 4D array with the same spatial dimensions as the input fields,
with the last dimension representing the 3 vector components.
"""
function poynting_vector(EMfield::ElectromagneticField{T}) where {T}
    # Ensure we have 4D arrays (3D space + 3 components)
    ndims(EMfield.E) == 4 ||
        throw(ArgumentError("poynting_vector requires 4D field arrays (3D space + 3 components)"))
    # Get field components
    Ex, Ey, Ez = EMfield.E[:, :, :, 1], EMfield.E[:, :, :, 2], EMfield.E[:, :, :, 3]
    Hx, Hy, Hz = EMfield.H[:, :, :, 1], EMfield.H[:, :, :, 2], EMfield.H[:, :, :, 3]

    # Calculate S = E × H* (complex conjugate of H)
    Sx = Ey .* conj.(Hz) .- Ez .* conj.(Hy)
    Sy = Ez .* conj.(Hx) .- Ex .* conj.(Hz)
    Sz = Ex .* conj.(Hy) .- Ey .* conj.(Hx)

    # Combine into 4D array
    S = similar(EMfield.E)
    S[:, :, :, 1] = Sx
    S[:, :, :, 2] = Sy
    S[:, :, :, 3] = Sz

    return S
end

"""
    field_intensity(EMfield::ElectromagneticField; component=:total) -> AbstractArray{T}

Calculate the field intensity |E|² or individual component intensities.

# Arguments
- `EMfield::ElectromagneticField`: Input electromagnetic fields
- `component::Symbol`: Component to calculate (:total, :x, :y, :z)

# Returns
- `intensity::AbstractArray{T}`: Field intensity distribution
"""
function field_intensity(EMfield::ElectromagneticField{T}; component::Symbol = :total) where {T}
    # Ensure we have 4D arrays (3D space + 3 components)
    ndims(EMfield.E) == 4 ||
        throw(ArgumentError("field_intensity requires 4D field arrays (3D space + 3 components)"))
    if component == :total
        return sum(abs2, EMfield.E; dims = 4)[:, :, :, 1]
    elseif component == :x
        return abs2.(EMfield.E[:, :, :, 1])
    elseif component == :y
        return abs2.(EMfield.E[:, :, :, 2])
    elseif component == :z
        return abs2.(EMfield.E[:, :, :, 3])
    else
        throw(ArgumentError("component must be :total, :x, :y, or :z"))
    end
end

"""
    extract_plane(EMfield::ElectromagneticField, plane_axis::Int, plane_index::Int) -> ElectromagneticField

Extract a 2D plane from the 3D electromagnetic field.

# Arguments
- `EMfield::ElectromagneticField`: Input 3D fields
- `plane_axis::Int`: Axis normal to the plane (1=X, 2=Y, 3=Z)
- `plane_index::Int`: Index along the plane axis

# Returns
- `plane_EMfield::ElectromagneticField`: 2D electromagnetic field on the specified plane
"""
function extract_plane(
        EMfield::ElectromagneticField{T},
        plane_axis::Int,
        plane_index::Int
) where {T}
    # Ensure we have 4D arrays (3D space + 3 components)
    ndims(EMfield.E) == 4 ||
        throw(ArgumentError("extract_plane requires 4D field arrays (3D space + 3 components)"))
    1 ≤ plane_axis ≤ 3 || throw(ArgumentError("plane_axis must be 1, 2, or 3"))
    1 ≤ plane_index ≤ EMfield.grid_size[plane_axis] ||
        throw(ArgumentError("plane_index out of bounds"))

    # Extract plane data
    if plane_axis == 1  # YZ plane
        E_plane = EMfield.E[plane_index:plane_index, :, :, :]
        H_plane = EMfield.H[plane_index:plane_index, :, :, :]
        plane_grid_size = (1, EMfield.grid_size[2], EMfield.grid_size[3])
        plane_resolution = (
            EMfield.resolution[1], EMfield.resolution[2], EMfield.resolution[3])
    elseif plane_axis == 2  # XZ plane
        E_plane = EMfield.E[:, plane_index:plane_index, :, :]
        H_plane = EMfield.H[:, plane_index:plane_index, :, :]
        plane_grid_size = (EMfield.grid_size[1], 1, EMfield.grid_size[3])
        plane_resolution = (
            EMfield.resolution[1], EMfield.resolution[2], EMfield.resolution[3])
    else  # XY plane
        E_plane = EMfield.E[:, :, plane_index:plane_index, :]
        H_plane = EMfield.H[:, :, plane_index:plane_index, :]
        plane_grid_size = (EMfield.grid_size[1], EMfield.grid_size[2], 1)
        plane_resolution = (
            EMfield.resolution[1], EMfield.resolution[2], EMfield.resolution[3])
    end

    return ElectromagneticField(
        E_plane, H_plane, plane_grid_size, plane_resolution, EMfield.wavelength)
end

"""
    +(field1::ElectromagneticField, field2::ElectromagneticField) -> ElectromagneticField

Add two electromagnetic fields element-wise.

Both fields must have the same grid size, resolution, and wavelength.
The resulting field contains the element-wise sum of E and H components.

# Example
```julia
total_field = scattered_field + incident_field
```
"""
function Base.:+(field1::ElectromagneticField{T}, field2::ElectromagneticField{T}) where {T}
    # Validate field compatibility
    field1.grid_size == field2.grid_size ||
        throw(ArgumentError("Fields must have same grid_size: $(field1.grid_size) vs $(field2.grid_size)"))
    field1.resolution == field2.resolution ||
        throw(ArgumentError("Fields must have same resolution"))
    isapprox(field1.wavelength, field2.wavelength, rtol=T(1e-10)) ||
        throw(ArgumentError("Fields must have same wavelength"))

    # Add E and H fields element-wise
    E_sum = field1.E .+ field2.E
    H_sum = field1.H .+ field2.H

    return ElectromagneticField(E_sum, H_sum, field1.grid_size, field1.resolution, field1.wavelength)
end

"""
    crop_to_ROI(field::ElectromagneticField, solver) -> ElectromagneticField

Crop electromagnetic field to the region of interest (ROI) defined by solver's boundary conditions.

This removes the padding that was added for boundary condition handling,
returning only the physical simulation domain.

# Arguments
- `field::ElectromagneticField{T}`: Padded electromagnetic field
- `solver`: Solver object containing ROI bounds

# Returns
- `ElectromagneticField{T}`: Cropped field on original grid

# Example
```julia
cropped_field = crop_to_ROI(padded_field, solver)
```
"""
function crop_to_ROI(field::ElectromagneticField{T}, solver) where {T}
    ROI = solver.ROI

    # Crop E and H arrays to ROI
    E_cropped = field.E[ROI[1]:ROI[2], ROI[3]:ROI[4], ROI[5]:ROI[6], :]
    H_cropped = field.H[ROI[1]:ROI[2], ROI[3]:ROI[4], ROI[5]:ROI[6], :]

    # Return cropped field with original grid size
    return ElectromagneticField(
        E_cropped,
        H_cropped,
        solver.grid_size,
        solver.resolution,
        field.wavelength
    )
end

"""
    to_host(field::ElectromagneticField) -> ElectromagneticField

Transfer electromagnetic field arrays from device (GPU) to host (CPU) memory.

This function handles data transfer for GPU-accelerated computations,
converting device arrays to standard Julia Arrays.

# Arguments
- `field::ElectromagneticField{T}`: Field with potentially device-resident arrays

# Returns
- `ElectromagneticField{T}`: Field with host-resident arrays

# Example
```julia
host_field = to_host(gpu_field)
```
"""
function to_host(field::ElectromagneticField{T}) where {T}
    # Transfer E and H to host
    E_host = Array{eltype(field.E)}(undef, size(field.E))
    H_host = Array{eltype(field.H)}(undef, size(field.H))

    copyto!(E_host, field.E)
    copyto!(H_host, field.H)

    return ElectromagneticField(E_host, H_host, field.grid_size, field.resolution, field.wavelength)
end

"""
    show(io::IO, EMfield::ElectromagneticField)

Custom display for electromagnetic field objects.
"""
function Base.show(io::IO, EMfield::ElectromagneticField{T}) where {T}
    print(io, "ElectromagneticField{$T}")
    print(io, "\n  grid_size: $(EMfield.grid_size)")
    print(io, "\n  resolution: $(EMfield.resolution)")
    print(io, "\n  domain_size: $(domain_size(EMfield))")
    print(io, "\n  wavelength: $(EMfield.wavelength)")
    print(io, "\n  field_energy: $(field_energy(EMfield))")
end
