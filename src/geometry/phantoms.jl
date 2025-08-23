"""
Phantom generation utilities for electromagnetic simulations.

This module provides functions to generate common geometric objects and
test phantoms for electromagnetic simulations, maintaining compatibility
with the legacy jl-ConvergentBornSolver API while adding type safety.
"""

"""
    phantom_bead(grid_size::NTuple{3, Int}, 
                 permittivity_profile::AbstractVector{<:Number}, 
                 pixel_radius::Real;
                 num_bead::Int = 1,
                 bead_distance::Real = 2 * pixel_radius,
                 array_type::Type{<:AbstractArray} = Array) -> AbstractArray{Complex{T}, 3}

Generate a phantom containing spherical beads with specified permittivity.

This function creates a 3D permittivity distribution containing one or more
spherical scatterers (beads) with specified optical properties.

# Arguments
- `grid_size::NTuple{3, Int}`: Grid dimensions (Nx, Ny, Nz)
- `permittivity_profile::AbstractVector{<:Number}`: Permittivity values for each bead
- `pixel_radius::Real`: Bead radius in grid pixels
- `num_bead::Int = 1`: Number of beads to generate
- `bead_distance::Real = 2 * pixel_radius`: Distance between bead centers (for multiple beads)
- `array_type::Type{<:AbstractArray} = Array`: Array type for output (CPU/GPU flexibility)

# Returns
- `phantom::AbstractArray{Complex{T}, 3}`: 3D permittivity distribution

# Algorithm
- Beads are placed along the Z-axis with specified spacing
- Each bead gets a permittivity value from `permittivity_profile` (cycling if needed)
- Background permittivity is 1.0 (vacuum/air)
- Spherical geometry: distance from bead center ≤ radius

# Example
```julia
# Single high-index bead
phantom = phantom_bead(
    (128, 128, 64),           # Grid size
    [1.5^2],                  # n = 1.5 refractive index
    8.0,                      # 8-pixel radius
    num_bead = 1
)

# Multiple beads with different indices
phantom = phantom_bead(
    (256, 256, 128),
    [1.3^2, 1.6^2, 1.4^2],    # Different materials
    10.0,                     # 10-pixel radius
    num_bead = 5,             # 5 beads total
    bead_distance = 25.0      # 25-pixel spacing
)
```
"""
function phantom_bead(
    grid_size::NTuple{3, Int},
    permittivity_profile::AbstractVector{<:Number},
    pixel_radius::Real;
    num_bead::Int = 1,
    bead_distance::Real = 2 * pixel_radius,
    array_type::Type{<:AbstractArray} = Array
)
    # Validate inputs
    all(grid_size .> 0) || throw(ArgumentError("grid_size must have positive dimensions"))
    pixel_radius > 0 || throw(ArgumentError("pixel_radius must be positive"))
    num_bead > 0 || throw(ArgumentError("num_bead must be positive"))
    !isempty(permittivity_profile) || throw(ArgumentError("permittivity_profile cannot be empty"))
    
    # Determine output type
    T = real(eltype(promote_type(eltype(permittivity_profile), typeof(pixel_radius))))
    T <: AbstractFloat || (T = Float64)
    
    # Initialize phantom with background permittivity
    phantom = array_type(ones(Complex{T}, grid_size))
    
    # Calculate grid center
    center = ntuple(i -> (grid_size[i] + 1) / 2, 3)
    
    # Place beads along Z-axis
    for bead_idx in 1:num_bead
        # Bead center position
        z_offset = (bead_idx - (num_bead + 1) / 2) * bead_distance
        bead_center = (center[1], center[2], center[3] + z_offset)
        
        # Ensure bead is within grid bounds
        if bead_center[3] < 1 || bead_center[3] > grid_size[3]
            @warn "Bead $bead_idx is outside grid bounds, skipping"
            continue
        end
        
        # Get permittivity for this bead (cycle through profile)
        perm_idx = ((bead_idx - 1) % length(permittivity_profile)) + 1
        bead_permittivity = Complex{T}(permittivity_profile[perm_idx])
        
        # Fill spherical region
        _fill_sphere!(phantom, bead_center, pixel_radius, bead_permittivity)
    end
    
    return phantom
end

"""
    phantom_plate(grid_size::NTuple{3, Int},
                  permittivity_profile::AbstractVector{<:Number},
                  thickness_pixel::Real;
                  plate_normal::Int = 3,
                  array_type::Type{<:AbstractArray} = Array) -> AbstractArray{Complex{T}, 3}

Generate a phantom containing a flat plate with specified permittivity.

This function creates a 3D permittivity distribution containing a flat plate
(slab) of specified thickness and optical properties.

# Arguments
- `grid_size::NTuple{3, Int}`: Grid dimensions (Nx, Ny, Nz)
- `permittivity_profile::AbstractVector{<:Number}`: Permittivity values (typically length 1)
- `thickness_pixel::Real`: Plate thickness in grid pixels
- `plate_normal::Int = 3`: Normal direction (1=X, 2=Y, 3=Z)
- `array_type::Type{<:AbstractArray} = Array`: Array type for output

# Returns
- `phantom::AbstractArray{Complex{T}, 3}`: 3D permittivity distribution

# Algorithm
- Plate is centered in the grid along the normal direction
- Extends fully in the two perpendicular directions
- Background permittivity is 1.0

# Example
```julia
# Glass plate in Z direction
phantom = phantom_plate(
    (128, 128, 64),           # Grid size
    [1.5^2],                  # Glass permittivity  
    8.0                       # 8-pixel thickness
)

# Thin film in Y direction
phantom = phantom_plate(
    (256, 128, 256),
    [2.1^2],                  # High-index material
    4.0,                      # 4-pixel thickness
    plate_normal = 2          # Y-normal
)
```
"""
function phantom_plate(
    grid_size::NTuple{3, Int},
    permittivity_profile::AbstractVector{<:Number},
    thickness_pixel::Real;
    plate_normal::Int = 3,
    array_type::Type{<:AbstractArray} = Array
)
    # Validate inputs
    all(grid_size .> 0) || throw(ArgumentError("grid_size must have positive dimensions"))
    thickness_pixel > 0 || throw(ArgumentError("thickness_pixel must be positive"))
    1 ≤ plate_normal ≤ 3 || throw(ArgumentError("plate_normal must be 1, 2, or 3"))
    !isempty(permittivity_profile) || throw(ArgumentError("permittivity_profile cannot be empty"))
    
    # Determine output type
    T = real(eltype(promote_type(eltype(permittivity_profile), typeof(thickness_pixel))))
    T <: AbstractFloat || (T = Float64)
    
    # Initialize phantom with background permittivity
    phantom = array_type(ones(Complex{T}, grid_size))
    
    # Get plate permittivity (use first value from profile)
    plate_permittivity = Complex{T}(permittivity_profile[1])
    
    # Calculate plate bounds along normal direction
    center_pos = (grid_size[plate_normal] + 1) / 2
    half_thickness = thickness_pixel / 2
    start_pos = max(1, round(Int, center_pos - half_thickness))
    end_pos = min(grid_size[plate_normal], round(Int, center_pos + half_thickness))
    
    # Fill plate region
    if plate_normal == 1  # X-normal plate
        phantom[start_pos:end_pos, :, :] .= plate_permittivity
    elseif plate_normal == 2  # Y-normal plate
        phantom[:, start_pos:end_pos, :] .= plate_permittivity
    else  # Z-normal plate (default)
        phantom[:, :, start_pos:end_pos] .= plate_permittivity
    end
    
    return phantom
end

"""
    _fill_sphere!(phantom::AbstractArray{Complex{T}, 3}, 
                  center::NTuple{3, Real}, 
                  radius::Real, 
                  permittivity::Complex{T}) where T

Fill a spherical region in the phantom array with specified permittivity.

This is an internal helper function that efficiently fills spherical regions
by iterating only over the bounding box and checking distances.
"""
function _fill_sphere!(
    phantom::AbstractArray{Complex{T}, 3},
    center::NTuple{3, Real},
    radius::Real,
    permittivity::Complex{T}
) where T
    
    grid_size = size(phantom)
    radius_sq = radius^2
    
    # Calculate bounding box to minimize iterations
    bbox_min = ntuple(i -> max(1, round(Int, center[i] - radius)), 3)
    bbox_max = ntuple(i -> min(grid_size[i], round(Int, center[i] + radius)), 3)
    
    # Fill spherical region
    for k in bbox_min[3]:bbox_max[3]
        dz = k - center[3]
        dz_sq = dz * dz
        
        for j in bbox_min[2]:bbox_max[2]
            dy = j - center[2]
            dy_sq = dy * dy
            
            for i in bbox_min[1]:bbox_max[1]
                dx = i - center[1]
                dx_sq = dx * dx
                
                # Check if point is inside sphere
                if dx_sq + dy_sq + dz_sq ≤ radius_sq
                    phantom[i, j, k] = permittivity
                end
            end
        end
    end
    
    return nothing
end

"""
    phantom_cylinder(grid_size::NTuple{3, Int},
                     permittivity_profile::AbstractVector{<:Number},
                     radius_pixel::Real;
                     height_pixel::Real = 2 * radius_pixel,
                     axis::Int = 3,
                     array_type::Type{<:AbstractArray} = Array) -> AbstractArray{Complex{T}, 3}

Generate a phantom containing a cylindrical scatterer.

# Arguments
- `grid_size::NTuple{3, Int}`: Grid dimensions
- `permittivity_profile::AbstractVector{<:Number}`: Cylinder permittivity values (typically length 1)
- `radius_pixel::Real`: Cylinder radius in pixels
- `height_pixel::Real = 2 * radius_pixel`: Cylinder height in pixels  
- `axis::Int = 3`: Cylinder axis direction (1=X, 2=Y, 3=Z)
- `array_type::Type{<:AbstractArray} = Array`: Array type for output

# Returns
- `phantom::AbstractArray{Complex{T}, 3}`: 3D permittivity distribution
"""
function phantom_cylinder(
    grid_size::NTuple{3, Int},
    permittivity_profile::AbstractVector{<:Number},
    radius_pixel::Real;
    height_pixel::Real = 2 * radius_pixel,
    axis::Int = 3,
    array_type::Type{<:AbstractArray} = Array
)
    # Validate inputs
    all(grid_size .> 0) || throw(ArgumentError("grid_size must have positive dimensions"))
    radius_pixel > 0 || throw(ArgumentError("radius_pixel must be positive"))
    height_pixel > 0 || throw(ArgumentError("height_pixel must be positive"))
    1 ≤ axis ≤ 3 || throw(ArgumentError("axis must be 1, 2, or 3"))
    !isempty(permittivity_profile) || throw(ArgumentError("permittivity_profile cannot be empty"))
    
    # Determine output type
    T = real(eltype(promote_type(eltype(permittivity_profile), typeof(radius_pixel))))
    T <: AbstractFloat || (T = Float64)
    
    # Initialize phantom
    phantom = array_type(ones(Complex{T}, grid_size))
    
    # Get cylinder permittivity (use first value from profile)
    cyl_permittivity = Complex{T}(permittivity_profile[1])
    
    # Get grid center
    center = ntuple(i -> (grid_size[i] + 1) / 2, 3)
    
    # Calculate cylinder bounds
    half_height = height_pixel / 2
    axis_start = max(1, round(Int, center[axis] - half_height))
    axis_end = min(grid_size[axis], round(Int, center[axis] + half_height))
    
    # Fill cylindrical region
    _fill_cylinder!(phantom, center, radius_pixel, axis_start, axis_end, axis, cyl_permittivity)
    
    return phantom
end
function _fill_cylinder!(
    phantom::AbstractArray{Complex{T}, 3},
    center::NTuple{3, Real},
    radius::Real,
    axis_start::Int,
    axis_end::Int,
    axis::Int,
    permittivity::Complex{T}
) where T
    
    radius_sq = radius^2
    
    if axis == 1  # X-axis cylinder
        for i in axis_start:axis_end
            for j in 1:size(phantom, 2), k in 1:size(phantom, 3)
                dy = j - center[2]
                dz = k - center[3]
                if dy*dy + dz*dz ≤ radius_sq
                    phantom[i, j, k] = permittivity
                end
            end
        end
    elseif axis == 2  # Y-axis cylinder
        for j in axis_start:axis_end
            for i in 1:size(phantom, 1), k in 1:size(phantom, 3)
                dx = i - center[1]
                dz = k - center[3]
                if dx*dx + dz*dz ≤ radius_sq
                    phantom[i, j, k] = permittivity
                end
            end
        end
    else  # Z-axis cylinder (default)
        for k in axis_start:axis_end
            for i in 1:size(phantom, 1), j in 1:size(phantom, 2)
                dx = i - center[1]
                dy = j - center[2]
                if dx*dx + dy*dy ≤ radius_sq
                    phantom[i, j, k] = permittivity
                end
            end
        end
    end
    
    return nothing
end
