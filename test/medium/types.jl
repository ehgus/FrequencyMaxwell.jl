"""
Medium type tests.
Tests Medium construction, accessors, validation, and type inference.
"""

@testset "Medium Type System" begin
    @testset "Basic Medium Construction" begin
        # Create a simple medium with homogeneous permittivity
        gs = (32, 32, 16)
        perm_dist = ones(ComplexF64, gs)
        perm_bg = 1.33^2  # Water background

        medium = Medium(perm_dist, perm_bg)

        # Test field access
        @test permittivity(medium) === perm_dist
        @test permittivity_bg(medium) ≈ perm_bg
        @test grid_size(medium) == gs
    end

    @testset "Type Parameter Inference" begin
        # Test Float64 construction
        perm_dist_f64 = ones(ComplexF64, 16, 16, 8)
        medium_f64 = Medium(perm_dist_f64, 1.0)
        @test eltype(permittivity_bg(medium_f64)) == Float64

        # Test Float32 construction
        perm_dist_f32 = ones(ComplexF32, 16, 16, 8)
        medium_f32 = Medium(perm_dist_f32, 1.0f0)
        @test eltype(permittivity_bg(medium_f32)) == Float32

        # Test automatic type conversion
        medium_auto = Medium(perm_dist_f64, 1)  # Integer background
        @test eltype(permittivity_bg(medium_auto)) == Float64
    end

    @testset "Accessor Methods" begin
        # Create medium with inhomogeneous permittivity
        gs = (20, 20, 10)
        perm_dist = ones(ComplexF64, gs)
        perm_dist[10:12, 10:12, 5:6] .= 2.25 + 0.1im  # Add scatterer
        perm_bg_val = 1.77  # Background refractive index ~1.33

        medium = Medium(perm_dist, perm_bg_val)

        # Test permittivity accessor
        ε = permittivity(medium)
        @test ε === perm_dist
        @test size(ε) == gs
        @test ε[10, 10, 5] ≈ 2.25 + 0.1im

        # Test background permittivity accessor
        ε_bg = permittivity_bg(medium)
        @test ε_bg ≈ perm_bg_val
        @test ε_bg isa Float64

        # Test grid size accessor
        retrieved_gs = grid_size(medium)
        @test retrieved_gs == gs
        @test retrieved_gs isa NTuple{3, Int}
    end

    @testset "Input Validation" begin
        # Test negative background permittivity (should fail)
        perm_dist = ones(ComplexF64, 10, 10, 10)
        @test_throws ArgumentError Medium(perm_dist, -1.0)
        @test_throws ArgumentError Medium(perm_dist, 0.0)

        # Test positive background permittivity (should pass)
        @test Medium(perm_dist, 1e-6) isa Medium
        @test Medium(perm_dist, 1.0) isa Medium
        @test Medium(perm_dist, 100.0) isa Medium
    end

    @testset "Different Grid Sizes" begin
        # Test various grid configurations
        grid_configs = [
            (8, 8, 8),      # Cubic
            (64, 64, 32),   # Typical simulation
            (100, 50, 25),  # Asymmetric
            (1, 1, 100),    # Thin slab
        ]

        for gs_config in grid_configs
            perm_dist = ones(ComplexF64, gs_config)
            medium = Medium(perm_dist, 1.0)
            @test grid_size(medium) == gs_config
        end
    end

    @testset "Complex Permittivity Distributions" begin
        # Test with realistic permittivity values
        gs = (32, 32, 16)

        # Lossless dielectric (real permittivity)
        perm_lossless = fill(2.25 + 0.0im, gs)
        medium_lossless = Medium(perm_lossless, 1.0)
        @test permittivity(medium_lossless)[1, 1, 1] ≈ 2.25 + 0.0im

        # Lossy dielectric (complex permittivity)
        perm_lossy = fill(2.25 + 0.5im, gs)
        medium_lossy = Medium(perm_lossy, 1.0)
        @test real(permittivity(medium_lossy)[1, 1, 1]) ≈ 2.25
        @test imag(permittivity(medium_lossy)[1, 1, 1]) ≈ 0.5

        # High-contrast scatterer
        perm_contrast = ones(ComplexF64, gs)
        perm_contrast[16, 16, 8] = 10.0 + 0.0im  # High permittivity inclusion
        medium_contrast = Medium(perm_contrast, 1.0)
        @test maximum(abs.(permittivity(medium_contrast))) ≈ 10.0
    end

    @testset "Show Method" begin
        # Test that show method works without errors
        gs = (64, 64, 32)
        perm_dist = ones(ComplexF64, gs)
        perm_dist[30:35, 30:35, 15:17] .= 2.25 + 0.1im
        perm_bg_val = 1.77

        medium = Medium(perm_dist, perm_bg_val)

        # Capture show output
        io = IOBuffer()
        show(io, medium)
        output = String(take!(io))

        # Verify output contains expected information
        @test occursin("Medium{Float64}", output)
        @test occursin("Grid size: 64 × 64 × 32", output)
        @test occursin("Background", output)
        @test occursin("εᵣ", output)
        @test occursin("Permittivity range", output)
        @test occursin("Refractive index range", output)
    end

    @testset "Type Stability" begin
        # Test that accessor functions are type-stable
        perm_dist = ones(ComplexF64, 16, 16, 8)
        medium = Medium(perm_dist, 1.77)

        # These should all be type-stable
        @inferred permittivity(medium)
        @inferred permittivity_bg(medium)
        @inferred grid_size(medium)
    end

    @testset "Medium Immutability" begin
        # Verify that Medium is immutable but its array content can change
        perm_dist = ones(ComplexF64, 10, 10, 10)
        medium = Medium(perm_dist, 1.0)

        # Modifying the underlying array should affect the medium
        perm_dist[5, 5, 5] = 2.0 + 0.0im
        @test permittivity(medium)[5, 5, 5] ≈ 2.0 + 0.0im

        # But we cannot reassign fields (should not compile if uncommented)
        # medium.permittivity_bg = 2.0  # This would error
    end

    @testset "Edge Cases" begin
        # Very small permittivity values
        perm_tiny = fill(1e-10 + 0.0im, 8, 8, 8)
        medium_tiny = Medium(perm_tiny, 1e-10)
        @test permittivity_bg(medium_tiny) > 0

        # Very large permittivity values
        perm_large = fill(1e6 + 0.0im, 8, 8, 8)
        medium_large = Medium(perm_large, 1e6)
        @test permittivity_bg(medium_large) ≈ 1e6

        # Mixed precision (Float32 with Float64 background)
        perm_f32 = ones(ComplexF32, 8, 8, 8)
        medium_mixed = Medium(perm_f32, 1.77)  # Float64 input
        @test permittivity_bg(medium_mixed) isa Float32  # Should convert to Float32
    end
end
