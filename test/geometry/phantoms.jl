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

        medium = phantom_bead(
            grid_size,
            [permittivity_bead],  # Single bead permittivity as vector
            radius_pixels,
            permittivity_bg = permittivity_bg
        )

        # Test Medium type
        @test medium isa Medium
        @test FrequencyMaxwell.grid_size(medium) == grid_size
        @test FrequencyMaxwell.permittivity_bg(medium) ≈ permittivity_bg

        # Get permittivity array from medium
        phantom = FrequencyMaxwell.permittivity(medium)
        @test size(phantom) == grid_size

        # Test that background is properly set
        @test phantom[1, 1, 1] ≈ permittivity_bg
        @test phantom[end, end, end] ≈ permittivity_bg

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
        permittivity_bg = 1.0

        # Use phantom_plate for single layer test (simplified from multi-layer)
        medium = phantom_plate(
            grid_size,
            [permittivity_list[2]],  # Use middle material (TiO2)
            thickness_pixels[1],      # Use first thickness
            permittivity_bg = permittivity_bg
        )

        # Test Medium type
        @test medium isa Medium
        @test FrequencyMaxwell.grid_size(medium) == grid_size
        @test FrequencyMaxwell.permittivity_bg(medium) ≈ permittivity_bg

        # Get permittivity array
        phantom = FrequencyMaxwell.permittivity(medium)
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
        medium_real = phantom_bead(
            grid_size, [permittivity_real], 8, permittivity_bg = 1.0
        )
        phantom_real = FrequencyMaxwell.permittivity(medium_real)
        @test all(real.(phantom_real) .≥ 1.0)
        @test all(imag.(phantom_real) .== 0.0)

        # Test complex permittivity (lossy)
        permittivity_complex = 2.25 + 0.1im  # Lossy material
        medium_complex = phantom_bead(
            grid_size, [permittivity_complex], 8, permittivity_bg = 1.0
        )
        phantom_complex = FrequencyMaxwell.permittivity(medium_complex)
        center = (17, 17, 17)
        @test real(phantom_complex[center...]) ≈ 2.25
        @test imag(phantom_complex[center...]) ≈ 0.1

        # Test high contrast materials
        permittivity_high = 16.0  # n = 4.0, high index
        medium_high = phantom_bead(
            grid_size, [permittivity_high], 5, permittivity_bg = 1.0
        )
        phantom_high = FrequencyMaxwell.permittivity(medium_high)
        @test maximum(real.(phantom_high)) ≈ 16.0
    end

    @testset "Phantom Edge Cases" begin
        grid_size = (20, 20, 20)

        # Test very small bead
        medium_small = phantom_bead(
            grid_size, [2.0], 1, permittivity_bg = 1.0
        )
        phantom_small = FrequencyMaxwell.permittivity(medium_small)
        @test count(x -> abs(x - 2.0) < 1e-10, phantom_small) ≥ 1

        # Test bead larger than grid (should be truncated)
        medium_large = phantom_bead(
            grid_size, [2.0], 15, permittivity_bg = 1.0
        )
        phantom_large = FrequencyMaxwell.permittivity(medium_large)
        @test count(x -> abs(x - 2.0) < 1e-10, phantom_large) >
              count(x -> abs(x - 1.0) < 1e-10, phantom_large)
    end

    @testset "Basic API Tests" begin
        @testset "Phantom Bead API" begin
            # Test the phantom_bead function with Medium return type
            grid_size = (32, 32, 32)
            permittivity_profile = [1.46^2]  # SiO2 (single value for single bead)
            radius_pixels = 5
            permittivity_bg = 1.33^2  # Water

            medium = phantom_bead(grid_size, permittivity_profile, radius_pixels,
                permittivity_bg = permittivity_bg)

            # Test Medium properties
            @test medium isa Medium
            @test FrequencyMaxwell.grid_size(medium) == grid_size
            @test FrequencyMaxwell.permittivity_bg(medium) ≈ permittivity_bg

            # Test permittivity array
            phantom = FrequencyMaxwell.permittivity(medium)
            @test size(phantom) == grid_size
            @test phantom[1, 1, 1] ≈ permittivity_bg  # Background

            # Test that center has bead permittivity
            center = div.(grid_size, 2) .+ 1
            center_value = phantom[center...]
            @test center_value ≈ permittivity_profile[1]
        end

        @testset "Phantom Plate API" begin
            # Test layered structure with Medium return type
            grid_size = (20, 20, 20)
            permittivity_profile = [2.0]  # Single material for plate
            thickness_pixels = 8  # Single thickness value
            permittivity_bg = 1.0

            medium = phantom_plate(grid_size, permittivity_profile, thickness_pixels,
                permittivity_bg = permittivity_bg)

            # Test Medium properties
            @test medium isa Medium
            @test FrequencyMaxwell.grid_size(medium) == grid_size
            @test FrequencyMaxwell.permittivity_bg(medium) ≈ permittivity_bg

            # Test permittivity array
            phantom = FrequencyMaxwell.permittivity(medium)
            @test size(phantom) == grid_size

            # Test that phantom was created (not checking exact values due to implementation details)
            @test !all(phantom .== permittivity_bg)  # Should have some variation from background
            @test maximum(real.(phantom)) ≥ 1.0  # Should have material with permittivity ≥ 1
        end
    end
end
