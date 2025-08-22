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
        
        phantom = create_spherical_phantom(
            grid_size, 
            permittivity_bg, 
            permittivity_bead, 
            radius_pixels
        )
        
        @test size(phantom) == grid_size
        
        # Test that background is properly set
        @test phantom[1, 1, 1] ≈ permittivity_bg
        @test phantom[end, end, end] ≈ permittivity_bg
        
        # Test that center contains bead material
        center = (div(grid_size[1], 2) + 1, div(grid_size[2], 2) + 1, div(grid_size[3], 2) + 1)
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
        
        phantom = create_layered_phantom(
            grid_size,
            permittivity_list,
            thickness_pixels
        )
        
        @test size(phantom) == grid_size
        
        # Test layer boundaries
        @test phantom[1, 25, 25] ≈ permittivity_list[1]  # First layer
        @test phantom[10, 25, 25] ≈ permittivity_list[2]  # Second layer  
        @test phantom[20, 25, 25] ≈ permittivity_list[3]  # Third layer
        
        # Test that layers are uniform in x-y plane
        for z in [5, 10, 15]
            ref_value = phantom[z, 25, 25]
            @test all(phantom[z, :, :] .≈ ref_value)
        end
    end
    
    @testset "Phantom Material Properties" begin
        # Test various material configurations
        grid_size = (32, 32, 32)
        
        # Test real permittivity (lossless)
        permittivity_real = 2.25  # n = 1.5
        phantom_real = create_spherical_phantom(
            grid_size, 1.0, permittivity_real, 8
        )
        @test all(real.(phantom_real) .≥ 1.0)
        @test all(imag.(phantom_real) .== 0.0)
        
        # Test complex permittivity (lossy)
        permittivity_complex = 2.25 + 0.1im  # Lossy material
        phantom_complex = create_spherical_phantom(
            grid_size, 1.0, permittivity_complex, 8
        )
        center = (17, 17, 17)
        @test real(phantom_complex[center...]) ≈ 2.25
        @test imag(phantom_complex[center...]) ≈ 0.1
        
        # Test high contrast materials
        permittivity_high = 16.0  # n = 4.0, high index
        phantom_high = create_spherical_phantom(
            grid_size, 1.0, permittivity_high, 5
        )
        @test maximum(real.(phantom_high)) ≈ 16.0
    end
    
    @testset "Phantom Edge Cases" begin
        grid_size = (20, 20, 20)
        
        # Test very small bead
        phantom_small = create_spherical_phantom(
            grid_size, 1.0, 2.0, 1
        )
        @test count(x -> abs(x - 2.0) < 1e-10, phantom_small) ≥ 1
        
        # Test bead larger than grid (should be truncated)
        phantom_large = create_spherical_phantom(
            grid_size, 1.0, 2.0, 15
        )
        @test count(x -> abs(x - 2.0) < 1e-10, phantom_large) > count(x -> abs(x - 1.0) < 1e-10, phantom_large)
        
        # Test zero radius (no bead)
        phantom_zero = create_spherical_phantom(
            grid_size, 1.5, 2.0, 0
        )
        @test all(phantom_zero .≈ 1.5)
    end
end
