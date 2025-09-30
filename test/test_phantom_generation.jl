"""
Phantom generation tests - migrated from forward_SiO2_5um_in_water.jl example
Tests geometric phantom creation and material property assignment.
"""

@testset "Phantom Generation" begin
    @testset "Spherical Bead Phantom" begin
        # Test basic bead phantom creation - matches original example
        grid_size = (101, 101, 96)  # Smaller than original for faster testing
        permittivity_bg = 1.333^2   # Water
        permittivity_bead = 1.4607^2  # SiO2
        radius_pixels = 10  # 0.5 micron at 50nm resolution

        phantom = phantom_bead(
            grid_size,
            [permittivity_bead],  # Single bead permittivity as vector
            radius_pixels
        )

        @test size(phantom) == grid_size

        # Test that background is properly set (phantom_bead uses 1.0 as background)
        @test phantom[1, 1, 1] ≈ 1.0
        @test phantom[end, end, end] ≈ 1.0

        # Test that center contains bead material
        center = (
            div(grid_size[1], 2) + 1, div(grid_size[2], 2) + 1, div(grid_size[3], 2) + 1)
        @test phantom[center...] ≈ permittivity_bead

        # Test phantom volume conservation
        bead_voxels = count(x -> abs(x - permittivity_bead) < 1e-10, phantom)
        expected_volume = (4/3) * π * radius_pixels^3
        @test bead_voxels / expected_volume ≈ 1.0 atol=0.3  # Allow for discretization error
    end

    @testset "Layered Phantom" begin
        # Test multi-layer phantom creation - based on grating example
        grid_size = (21, 50, 50)  # Smaller for testing
        permittivity_list = [1.4338^2, (2.9734 + 0.0467im)^2, 1.6380^2]  # PDMS, TiO2, SU8
        thickness_pixels = [4, 15]  # Layer thicknesses

        # Use phantom_plate for single layer test (simplified from multi-layer)
        phantom = phantom_plate(
            grid_size,
            [permittivity_list[2]],  # Use middle material (TiO2)
            thickness_pixels[1]      # Use first thickness
        )

        @test size(phantom) == grid_size

        # Test that plate is created (simplified test for single layer)
        @test size(phantom) == grid_size

        # Test that phantom has the expected material in the center region
        center_z = div(grid_size[3], 2) + 1
        @test phantom[10, 25, center_z] ≈ permittivity_list[2]  # Should have TiO2 material

        # Test that plate is uniform in x-y plane at center
        ref_value = phantom[10, 25, center_z]
        @test all(phantom[10, :, center_z] .≈ ref_value)
    end

    @testset "Phantom Material Properties" begin
        # Test various material configurations
        grid_size = (32, 32, 32)

        # Test real permittivity (lossless)
        permittivity_real = 2.25  # n = 1.5
        phantom_real = phantom_bead(
            grid_size, [permittivity_real], 8
        )
        @test all(real.(phantom_real) .≥ 1.0)
        @test all(imag.(phantom_real) .== 0.0)

        # Test complex permittivity (lossy)
        permittivity_complex = 2.25 + 0.1im  # Lossy material
        phantom_complex = phantom_bead(
            grid_size, [permittivity_complex], 8
        )
        center = (17, 17, 17)
        @test real(phantom_complex[center...]) ≈ 2.25
        @test imag(phantom_complex[center...]) ≈ 0.1

        # Test high contrast materials
        permittivity_high = 16.0  # n = 4.0, high index
        phantom_high = phantom_bead(
            grid_size, [permittivity_high], 5
        )
        @test maximum(real.(phantom_high)) ≈ 16.0
    end

    @testset "Phantom Edge Cases" begin
        grid_size = (20, 20, 20)

        # Test very small bead
        phantom_small = phantom_bead(
            grid_size, [2.0], 1
        )
        @test count(x -> abs(x - 2.0) < 1e-10, phantom_small) ≥ 1

        # Test bead larger than grid (should be truncated)
        phantom_large = phantom_bead(
            grid_size, [2.0], 15
        )
        @test count(x -> abs(x - 2.0) < 1e-10, phantom_large) >
              count(x -> abs(x - 1.0) < 1e-10, phantom_large)
    end
end
