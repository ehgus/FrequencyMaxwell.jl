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
struct Curl{T <: Real}
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
    conv!(curl::Curl, field::AbstractArray) -> AbstractArray

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
function conv!(curl::Curl, field::AbstractArray)
    # Transform to frequency domain
    fft!(field, 1:3)

    # Apply curl function in frequency domain
    multiply_curl!(curl, field)

    # Transform back to spatial domain
    ifft!(field, 1:3)

    return field
end

"""
GPU kernel for curl operation in frequency domain.

Computes the curl operation: (∇×) = (ik ×)
"""
@kernel function multiply_curl_kernel!(field, freq_resolution)
    I = @index(Global, Cartesian)
    i, j, k = I[1], I[2], I[3]
    max_i, max_j, max_k = @ndrange()

    # Ensure bounds (safety check)
    if i <= max_i && j <= max_j && k <= max_k
        src1 = field[i, j, k, 1] # Ex
        src2 = field[i, j, k, 2] # Ey
        src3 = field[i, j, k, 3] # Ez

        idx1 = i <= div(max_i, 2) ? i - 1 : i - max_i - 1
        idx2 = j <= div(max_j, 2) ? j - 1 : j - max_j - 1
        idx3 = k <= div(max_k, 2) ? k - 1 : k - max_k - 1

        ikx = 1im * freq_resolution[1] * idx1
        iky = 1im * freq_resolution[2] * idx2
        ikz = 1im * freq_resolution[3] * idx3

        dst1 = iky * src3 - ikz * src2 # (∇ × E)_x
        dst2 = ikz * src1 - ikx * src3 # (∇ × E)_y
        dst3 = ikx * src2 - iky * src1 # (∇ × E)_z
        field[i, j, k, 1] = dst1
        field[i, j, k, 2] = dst2
        field[i, j, k, 3] = dst3
    end
end

function multiply_curl!(curl::Curl, field::AbstractArray)
    freq_resolution = curl.freq_res

    # Get appropriate backend for hardware
    backend = get_backend(field)

    kernel! = multiply_curl_kernel!(backend, 64)
    kernel!(field, freq_resolution, ndrange = size(field)[1:3])
end

function conv(curl::Curl, field::AbstractArray)
    backend = get_backend(field)
    result = KernelAbstractions.zeros(backend, eltype(field), size(field))
    result .= field
    conv!(curl, result)
end
