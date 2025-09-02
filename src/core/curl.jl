"""
Curl operator implementation for electromagnetic field computation.

This module provides efficient curl computation in Fourier space
for electromagnetic field analysis.
"""

using FFTW

"""
    Curl{T}

Curl operator for electromagnetic fields.

# Fields
- `array_type::Type{<:AbstractArray}`: Array type for hardware flexibility
- `precision::Type{<:Real}`: Floating-point precision type
- `freq_res::NTuple{3, T}`: Frequency resolution for each dimension
"""
struct Curl{T<:Real}
    array_type::Type{<:AbstractArray}
    precision::Type{<:Real}
    freq_res::NTuple{3, T}
end

"""
    Curl(array_type, precision, arr_size, resolution)

Construct a curl operator.

# Arguments
- `array_type`: Array type (Array, CuArray, etc.)
- `precision`: Floating-point precision (Float32, Float64)
- `arr_size`: Domain size (nx, ny, nz)
- `resolution`: Spatial resolution (dx, dy, dz)
"""
function Curl(array_type::Type, precision::Type, arr_size::NTuple{3, <:Integer}, 
              resolution::NTuple{3, <:Real})
    T = precision <: Real ? precision : Float64
    # Calculate frequency resolution
    freq_res = 2π ./ (arr_size .* resolution)
    return Curl{T}(array_type, precision, T.(freq_res))
end

"""
    conv(curl::Curl, field::AbstractArray) -> AbstractArray

Compute curl of electromagnetic field using FFT.

The curl operation is performed in Fourier space:
∇ × E → ik × Ê

where k is the wave vector and Ê is the Fourier transform of E.

# Arguments
- `curl::Curl`: Curl operator
- `field::AbstractArray`: Input electromagnetic field (nx, ny, nz, 3)

# Returns
- `AbstractArray`: Curl of the field with same dimensions
"""
function conv(curl::Curl, field::AbstractArray)
    # Initialize result array
    result = similar(field)
    fill!(result, zero(eltype(result)))

    fourier_resolution = curl.freq_res
    
    # Precompute frequency coordinates for each dimension
    fourier_coord = Vector{Any}()
    
    for axis = 1:3
        # Create frequency coordinates: [0, 1, 2, ..., N/2-1, -N/2, -N/2+1, ..., -1]
        N = size(result, axis)
        coor_axis = vcat(0:div(N,2)-1, -div(N,2):-1)
        if isodd(N)
            coor_axis = vcat(0:div(N,2), -div(N,2):-1)
        end
        
        # Scale by frequency resolution and multiply by i for derivative
        coor_axis = coor_axis .* fourier_resolution[axis]
        
        # Reshape for broadcasting: put the axis dimension in correct position
        dims = ones(Int, 3)
        dims[axis] = N
        coor_axis_shaped = reshape(coor_axis, Tuple(dims))
        
        # Convert to complex and multiply by i for derivative  
        coord_complex = Complex{curl.precision}(1im) * coor_axis_shaped
        push!(fourier_coord, coord_complex)
    end
    
    # Compute curl in Fourier space: ∇ × E = ik × Ê
    # For each component: (∇ × E)_i = ε_ijk * ik_j * Ê_k
    
    for axis in 1:3
        # Get cyclic permutation for cross product
        # axis=1: (2,3) -> y,z; axis=2: (3,1) -> z,x; axis=3: (1,2) -> x,y
        axes_ordered = [mod1(axis, 3), mod1(axis+1, 3), mod1(axis+2, 3)]
        
        # Take FFT of the field component we're operating on
        field_component_1 = fft(field[:, :, :, axes_ordered[3]])  # Third component
        field_component_2 = fft(field[:, :, :, axes_ordered[2]])  # Second component
        
        # Apply curl: (∇ × E)_i = ik_j * E_k - ik_k * E_j
        result[:, :, :, axes_ordered[1]] .+= fourier_coord[axes_ordered[2]] .* field_component_1
        result[:, :, :, axes_ordered[1]] .-= fourier_coord[axes_ordered[3]] .* field_component_2
    end
    
    # Transform back to spatial domain
    ifft!(result, 1:3)
    
    return result
end