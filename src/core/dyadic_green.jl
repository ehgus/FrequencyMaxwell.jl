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
- `shifts::NTuple{3, T}`: Subpixel shift for boundary condition handling
- `use_averaging::Bool`: Whether to use (G + G_flip)/2 averaging for boundary conditions

# Boundary Condition Handling

The `use_averaging` field controls CBS algorithm behavior:
- `false`: Periodic boundaries - G and G_flip identical (both use zero shift)
- `true`: Non-periodic boundaries - Average G(+shift) and G(-shift)

When `use_averaging = true`, the `conv_with_averaging!` method internally:
1. Applies G with positive subpixel shift
2. Applies G with negative subpixel shift (flip)
3. Averages the results: (G + G_flip)/2

This eliminates the need for manual flip management in the solver.

# Constructor
```julia
bc_periodic = PeriodicBoundaryCondition()
bc_absorbing = AbsorbingBoundaryCondition(thickness=3e-6, attenuation_thickness=2e-6)

green = DyadicGreen(Array, (128, 128, 64), k², (50e-9, 50e-9, 50e-9),
                   (bc_periodic, bc_periodic, bc_absorbing))
```
"""
struct DyadicGreen{T <: Real}
    flip_array::AbstractArray
    k_square::Complex{T}
    freq_res::NTuple{3, T}
    shifts::NTuple{3, T}
    use_averaging::Bool

    """
    DyadicGreen(array_type, arr_size, k_square, resolution, boundary_conditions)

    Construct dyadic Green's function operator from boundary conditions.

    This constructor automatically extracts subpixel shifts and averaging requirements
    from the provided boundary conditions.

    # Arguments
    - `array_type::Type`: Array type for GPU/CPU backend (e.g., Array, CuArray)
    - `k_square::Number`: Background wavenumber squared (k₀²)
    - `arr_size::NTuple{3, Integer}`: Padded grid dimensions (nx, ny, nz)
    - `resolution::NTuple{3, Real}`: Spatial resolution (dx, dy, dz) in meters
    - `boundary_conditions::NTuple{3, AbstractBoundaryCondition}`: Boundary conditions per dimension
    """
    function DyadicGreen(
            backend::Backend, arr_size::NTuple{3, <:Integer},
            k_square::Number, resolution::NTuple{3, <:Real},
            boundary_conditions::Tuple{Vararg{<:AbstractBoundaryCondition, 3}})
        # Extract subpixel shifts from boundary conditions
        shifts = ntuple(3) do i
            subpixel_shift(boundary_conditions[i])
        end

        # Check if any boundary condition requires averaging
        use_averaging = any(requires_averaging, boundary_conditions)
        T = real(typeof(k_square))

        # Allocate tmp_array when using averaging
        if use_averaging
            field_size = (arr_size[1:3]..., 3)
            flip_array = KernelAbstractions.zeros(backend, Complex{T}, field_size)
        else
            flip_array = []
        end

        # Calculate frequency resolution: Δk = 2π/(N*Δx)
        freq_res = 2π ./ (arr_size .* resolution)
        return new{T}(flip_array, Complex{T}(k_square), T.(freq_res),
            T.(shifts), use_averaging)
    end
end

function DyadicGreen(arr::AbstractArray, k_square, resolution, boundary_conditions)
    backend = get_backend(arr)
    arr_size = size(arr)
    return DyadicGreen(backend, arr_size, k_square, resolution, boundary_conditions)
end

"""
    conv!(green::DyadicGreen, src::AbstractArray)

In-place convolution with dyadic Green's function including automatic flip averaging.

Convolution with dyadic Green's function: src = G ⊛ src

If `green.use_averaging == true`, automatically computes (G + G_flip)/2 by averaging
convolutions with positive and negative subpixel shifts. This is required for non-periodic
boundary conditions in the CBS algorithm.

The convolution is performed in Fourier space for efficiency:
1. Apply phase ramp (for subpixel shifts)
2. Forward FFT
3. Multiply by Green's function in frequency domain
4. Inverse FFT
5. Apply inverse phase ramp
6. If averaging: repeat with negated shift and average results

# Arguments
- `green::DyadicGreen`: Green's function operator
- `src::AbstractArray`: Source array (modified in-place)

# Returns
- Modified `src` array
"""
function conv!(green::DyadicGreen, src::AbstractArray)
    if !green.use_averaging
        # Simple convolution without averaging (periodic boundaries)
        _apply_single_convolution!(green, src, green.shifts)
        return src
    else
        # Convolution with averaging: (G + G_flip)/2 (non-periodic boundaries)

        # Save original for flip convolution
        flip_src = green.flip_array
        flip_src .= src

        # Apply G with positive shift
        _apply_single_convolution!(green, src, green.shifts)

        # Apply G with negative shift (flip)
        negated_shift = .- green.shifts
        _apply_single_convolution!(green, flip_src, negated_shift)

        # Average: (G + G_flip)/2
        src .+= flip_src
        src ./= 2

        return src
    end
end

"""
    _apply_single_convolution!(green::DyadicGreen, src::AbstractArray, shift::NTuple{3})

Internal method: Apply a single Green's function convolution with specified subpixel shift.

This performs the core convolution operation in Fourier space without averaging.
"""
function _apply_single_convolution!(green::DyadicGreen, src::AbstractArray,
        shifts::NTuple{3, <:Real})
    # Phase ramp for subpixel boundary conditions
    _multiply_phase_ramp!(src, shifts, false)

    # Transform to frequency domain
    fft!(src, 1:3)

    # Apply Dyadic Green's function in frequency domain
    _multiply_Green!(green, src, shifts)

    # Transform back to spatial domain
    ifft!(src, 1:3)

    # Inverse phase ramp
    _multiply_phase_ramp!(src, shifts, true)

    return src
end

"""
GPU kernel for dyadic Green's function multiplication in frequency domain.

Computes the electromagnetic dyadic Green's function:
Ĝ(k) = (1/k²)[I - kk/k²] = (1/(k² - k²ₘ))[I - kk/k²]

where k² = k²ₓ + k²ᵧ + k²ᵤ and k²ₘ is the background wave number.
"""
@kernel function multiply_Green_kernel!(freq_res, shifts, k_square, src)
    I = @index(Global, Cartesian)
    i, j, k = I[1], I[2], I[3]
    max_i, max_j, max_k = @ndrange()

    # Ensure bounds (safety check)
    if i <= max_i && j <= max_j && k <= max_k
        # Extract field components
        src1 = src[i, j, k, 1]  # Ex
        src2 = src[i, j, k, 2]  # Ey  
        src3 = src[i, j, k, 3]  # Ez

        # Compute k-vector components with proper frequency centering
        # k = 2π * freq_index / domain_size, with DC at center
        kx_val = freq_res[1] *
                 (rem(i + floor(max_i/2) - 1, max_i) - floor(max_i/2) + shifts[1])
        ky_val = freq_res[2] *
                 (rem(j + floor(max_j/2) - 1, max_j) - floor(max_j/2) + shifts[2])
        kz_val = freq_res[3] *
                 (rem(k + floor(max_k/2) - 1, max_k) - floor(max_k/2) + shifts[3])

        # Compute k·E (divergence in frequency domain)
        k_dot_E = (kx_val * src1 + ky_val * src2 + kz_val * src3) / k_square

        # Compute denominator: k² - k²_medium
        k_square_diff = abs(kx_val^2 + ky_val^2 + kz_val^2) - k_square

        # Apply dyadic Green's function: Ĝ·E = (1/(k² - k²ₘ))[E - k(k·E)/k²]
        src[i, j, k, 1] = (src1 - kx_val * k_dot_E) / k_square_diff
        src[i, j, k, 2] = (src2 - ky_val * k_dot_E) / k_square_diff
        src[i, j, k, 3] = (src3 - kz_val * k_dot_E) / k_square_diff
    end
end

"""
Internal: Launch dyadic Green's function multiplication kernel with specified shift.
"""
function _multiply_Green!(green::DyadicGreen, src::AbstractArray, shifts::NTuple{3, <:Real})
    freq_res = green.freq_res
    k_square = green.k_square

    # Get appropriate backend for hardware
    backend = get_backend(src)

    # Launch kernel with block size 64 (good for both CPU and GPU)
    kernel! = multiply_Green_kernel!(backend, 64)
    kernel!(freq_res, shifts, k_square, src, ndrange = size(src)[1:3])
end

"""
GPU kernel for phase ramp multiplication.

Applies exp(±i·k·r) phase ramps for subpixel boundary condition handling.
"""
@kernel function multiply_phase_ramp_kernel!(arr::AbstractArray, shifts, do_reverse::Bool)
    I = @index(Global, Cartesian)
    i, j, k = I[1], I[2], I[3]
    max_i, max_j, max_k = @ndrange()

    if i <= max_i && j <= max_j && k <= max_k
        # Compute phase for subpixel shift
        phase = (i - 1) * shifts[1] / max_i +
                (j - 1) * shifts[2] / max_j +
                (k - 1) * shifts[3] / max_k

        # Apply forward or reverse phase ramp
        complex_phase = if do_reverse
            exp(2π * im * phase)
        else
            exp(-2π * im * phase)
        end

        # Apply phase to all field components
        for axis in 1:3
            arr[i, j, k, axis] *= complex_phase
        end
    end
end

"""
Internal: Apply phase ramp for subpixel boundary condition handling with specified shift.
"""
function _multiply_phase_ramp!(arr::AbstractArray, shifts::NTuple{3, <:Real}, do_reverse::Bool)
    # Get backend and launch kernel
    backend = get_backend(arr)
    kernel! = multiply_phase_ramp_kernel!(backend, 64)
    kernel!(arr, shifts, do_reverse, ndrange = size(arr)[1:3])
end
