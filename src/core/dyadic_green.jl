"""
Dyadic Green's function implementation for electromagnetic scattering.

This module provides GPU-accelerated dyadic Green's function computation
using KernelAbstractions.jl for hardware portability.
"""

using KernelAbstractions
using FFTW

"""
    DyadicGreen{T}

Dyadic Green's function operator for electromagnetic field computation.

# Fields
- `array_type::Type{<:AbstractArray}`: Array type for hardware flexibility (Array, CuArray, etc.)
- `k_square::Complex{T}`: Squared wave number in the medium
- `freq_res::NTuple{3, T}`: Frequency resolution for each spatial dimension
- `subpixel_shift::NTuple{3, T}`: Subpixel shift for boundary condition handling
"""
struct DyadicGreen{T<:Real}
    array_type::Type{<:AbstractArray}
    k_square::Complex{T}
    freq_res::NTuple{3, T}
    subpixel_shift::NTuple{3, T}
end

"""
    DyadicGreen(array_type, k_square, arr_size, resolution, subpixel_shift)

Construct a DyadicGreen operator.

# Arguments
- `array_type`: Array type for computation (Array, CuArray, etc.)
- `k_square`: Squared wave number k² = (ω/c)²ε_bg
- `arr_size`: Size of computational domain (nx, ny, nz)
- `resolution`: Spatial resolution (dx, dy, dz)
- `subpixel_shift`: Subpixel shift for boundary conditions
"""
function DyadicGreen(array_type::Type, k_square::Number, arr_size::NTuple{3, <:Integer}, 
                    resolution::NTuple{3, <:Real}, subpixel_shift::NTuple{3, <:Real})
    T = real(typeof(k_square))
    # Calculate frequency resolution: Δk = 2π/(N*Δx)
    freq_res = 2π ./ (arr_size .* resolution)
    return DyadicGreen{T}(array_type, Complex{T}(k_square), T.(freq_res), T.(subpixel_shift))
end

"""
    conv!(green::DyadicGreen, src::AbstractArray)

In-place convolution with dyadic Green's function.
"""
function conv!(green::DyadicGreen, src::AbstractArray)
    conv!(green, src, src)
end

"""
    conv!(green::DyadicGreen, src::AbstractArray, dst::AbstractArray)

Convolution with dyadic Green's function: dst = G ⊛ src

The convolution is performed in Fourier space for efficiency:
1. Apply phase ramp (for subpixel shifts)
2. Forward FFT
3. Multiply by Green's function in frequency domain
4. Inverse FFT
5. Apply inverse phase ramp
"""
function conv!(green::DyadicGreen, src::AbstractArray, dst::AbstractArray)
    # Phase ramp for subpixel boundary conditions
    multiply_phase_ramp!(green, src, do_reverse = false)
    
    # Transform to frequency domain
    fft!(src, 1:3)
    
    # Apply Dyadic Green's function in frequency domain
    multiply_Green!(green, dst, src)
    
    # Transform back to spatial domain
    ifft!(dst, 1:3)
    
    # Inverse phase ramp
    multiply_phase_ramp!(green, dst, do_reverse = true)
    
    return dst
end

"""
GPU kernel for dyadic Green's function multiplication in frequency domain.

Computes the electromagnetic dyadic Green's function:
Ĝ(k) = (1/k²)[I - kk/k²] = (1/(k² - k²ₘ))[I - kk/k²]

where k² = k²ₓ + k²ᵧ + k²ᵤ and k²ₘ is the background wave number.
"""
@kernel function multiply_Green_kernel!(freq_res, subpixel_shift, k_square, src, dst)
    I = @index(Global, Cartesian)
    i, j, k = I[1], I[2], I[3]
    max_i, max_j, max_k = @ndrange()

    # Ensure bounds (safety check)
    if i <= max_i && j <= max_j && k <= max_k
        # Extract field components
        src1 = src[i,j,k,1]  # Ex
        src2 = src[i,j,k,2]  # Ey  
        src3 = src[i,j,k,3]  # Ez

        # Compute k-vector components with proper frequency centering
        # k = 2π * freq_index / domain_size, with DC at center
        kx_val = freq_res[1] * (rem(i + floor(max_i/2) - 1, max_i) - floor(max_i/2) + subpixel_shift[1])
        ky_val = freq_res[2] * (rem(j + floor(max_j/2) - 1, max_j) - floor(max_j/2) + subpixel_shift[2])
        kz_val = freq_res[3] * (rem(k + floor(max_k/2) - 1, max_k) - floor(max_k/2) + subpixel_shift[3])

        # Compute k·E (divergence in frequency domain)
        k_dot_E = (kx_val * src1 + ky_val * src2 + kz_val * src3) / k_square

        # Compute denominator: k² - k²_medium
        k_square_diff = abs(kx_val^2 + ky_val^2 + kz_val^2) - k_square

        # Apply dyadic Green's function: Ĝ·E = (1/(k² - k²ₘ))[E - k(k·E)/k²]
        dst[i,j,k,1] = (src1 - kx_val * k_dot_E) / k_square_diff
        dst[i,j,k,2] = (src2 - ky_val * k_dot_E) / k_square_diff
        dst[i,j,k,3] = (src3 - kz_val * k_dot_E) / k_square_diff
    end
end

"""
Launch dyadic Green's function multiplication kernel.
"""
function multiply_Green!(green::DyadicGreen, dst::AbstractArray, src::AbstractArray)
    freq_res = green.freq_res
    subpixel_shift = green.subpixel_shift
    k_square = green.k_square

    # Get appropriate backend for hardware
    backend = get_backend(src)
    
    # Launch kernel with block size 64 (good for both CPU and GPU)
    kernel! = multiply_Green_kernel!(backend, 64)
    kernel!(freq_res, subpixel_shift, k_square, src, dst, 
            ndrange = size(src)[1:3])
end

"""
GPU kernel for phase ramp multiplication.

Applies exp(±i·k·r) phase ramps for subpixel boundary condition handling.
"""
@kernel function multiply_phase_ramp_kernel!(arr::AbstractArray, subpixel_shift, do_reverse::Bool)
    I = @index(Global, Cartesian)
    i, j, k = I[1], I[2], I[3]
    max_i, max_j, max_k = @ndrange()

    if i <= max_i && j <= max_j && k <= max_k
        # Compute phase for subpixel shift
        phase = (i - 1) * subpixel_shift[1] / max_i +
                (j - 1) * subpixel_shift[2] / max_j +
                (k - 1) * subpixel_shift[3] / max_k

        # Apply forward or reverse phase ramp
        complex_phase = if do_reverse
            exp(2π * im * phase)
        else
            exp(-2π * im * phase)
        end

        # Apply phase to all field components
        for axis = 1:3
            arr[i,j,k,axis] *= complex_phase
        end
    end
end

"""
Apply phase ramp for subpixel boundary condition handling.
"""
function multiply_phase_ramp!(green::DyadicGreen, arr::AbstractArray; do_reverse = false)
    subpixel_shift = green.subpixel_shift
    
    # Get backend and launch kernel
    backend = get_backend(arr)
    kernel! = multiply_phase_ramp_kernel!(backend, 64)
    kernel!(arr, subpixel_shift, do_reverse, 
            ndrange = size(arr)[1:3])
end