"""
Boundary condition tests - comprehensive interface and functionality tests
Tests the boundary condition abstraction and its integration with solvers.
"""

@testset "Boundary Condition Interface" begin
    @testset "PeriodicBoundaryCondition Construction" begin
        # Test default Float64 construction
        bc_default = PeriodicBoundaryCondition()
        @test bc_default isa PeriodicBoundaryCondition{Float64}

        # Test explicit type construction
        bc_f32 = PeriodicBoundaryCondition(Float32)
        @test bc_f32 isa PeriodicBoundaryCondition{Float32}

        bc_f64 = PeriodicBoundaryCondition(Float64)
        @test bc_f64 isa PeriodicBoundaryCondition{Float64}
    end

    @testset "AbsorbingBoundaryCondition Construction" begin
        # Test basic construction with defaults
        bc = AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.5e-6
        )

        @test bc isa AbsorbingBoundaryCondition
        @test bc.thickness ≈ 2.0e-6
        @test bc.attenuation_thickness ≈ 1.5e-6
        @test bc.sharpness ≈ 0.9  # Default
        @test bc.profile == TanhProfile  # Default

        # Test with explicit parameters
        bc_explicit = AbsorbingBoundaryCondition(
            thickness = 3.0e-6,
            attenuation_thickness = 2.0e-6,
            sharpness = 0.85,
            profile = ExponentialProfile
        )

        @test bc_explicit.thickness ≈ 3.0e-6
        @test bc_explicit.sharpness ≈ 0.85
        @test bc_explicit.profile == ExponentialProfile
    end

    @testset "AbsorbingBoundaryCondition Validation" begin
        # Test negative thickness
        @test_throws ArgumentError AbsorbingBoundaryCondition(
            thickness = -1.0e-6,
            attenuation_thickness = 0.5e-6
        )

        # Test negative attenuation_thickness
        @test_throws ArgumentError AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = -0.5e-6
        )

        # Test attenuation_thickness > thickness
        @test_throws ArgumentError AbsorbingBoundaryCondition(
            thickness = 1.0e-6,
            attenuation_thickness = 2.0e-6
        )

        # Test sharpness < 0
        @test_throws ArgumentError AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.0e-6,
            sharpness = -0.1
        )

        # Test sharpness > 1
        @test_throws ArgumentError AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.0e-6,
            sharpness = 1.5
        )
    end

    @testset "padding_pixels Interface" begin
        # Test periodic boundary (no padding)
        bc_periodic = PeriodicBoundaryCondition()
        @test FrequencyMaxwell.padding_pixels(bc_periodic, 50e-9) == 0
        @test FrequencyMaxwell.padding_pixels(bc_periodic, 100e-9) == 0

        # Test absorbing boundary
        bc_absorbing = AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.5e-6
        )

        # Should return thickness / (resolution * 2)
        @test FrequencyMaxwell.padding_pixels(bc_absorbing, 50e-9) ==
              round(Int, 2.0e-6 / (50e-9 * 2))
        @test FrequencyMaxwell.padding_pixels(bc_absorbing, 100e-9) ==
              round(Int, 2.0e-6 / (100e-9 * 2))
    end

    @testset "subpixel_shift Interface" begin
        # Test periodic boundary (zero shift)
        bc_periodic = PeriodicBoundaryCondition()
        @test FrequencyMaxwell.subpixel_shift(bc_periodic) == 0.0

        # Test absorbing boundary (0.25 shift)
        bc_absorbing = AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.5e-6
        )
        @test FrequencyMaxwell.subpixel_shift(bc_absorbing) == 0.25
    end

    @testset "requires_padding Interface" begin
        # Test periodic boundary
        bc_periodic = PeriodicBoundaryCondition()
        @test FrequencyMaxwell.requires_padding(bc_periodic) == false

        # Test absorbing boundary
        bc_absorbing = AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.5e-6
        )
        @test FrequencyMaxwell.requires_padding(bc_absorbing) == true
    end

    @testset "requires_averaging Interface" begin
        # Test periodic boundary (no averaging needed)
        bc_periodic = PeriodicBoundaryCondition()
        @test FrequencyMaxwell.requires_averaging(bc_periodic) == false

        # Test absorbing boundary (averaging required)
        bc_absorbing = AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.5e-6
        )
        @test FrequencyMaxwell.requires_averaging(bc_absorbing) == true
    end
end

@testset "Attenuation Profiles" begin
    @testset "Profile Enumeration" begin
        # Test that all profile types are defined
        @test isdefined(FrequencyMaxwell, :TanhProfile)
        @test isdefined(FrequencyMaxwell, :ExponentialProfile)
        @test isdefined(FrequencyMaxwell, :TukeyProfile)
    end

    @testset "Profile Selection" begin
        profiles = [TanhProfile, ExponentialProfile, TukeyProfile]

        for profile in profiles
            bc = AbsorbingBoundaryCondition(
                thickness = 2.0e-6,
                attenuation_thickness = 1.5e-6,
                profile = profile
            )

            @test bc.profile == profile
        end
    end

    @testset "Profile Default" begin
        # Test that TanhProfile is the default
        bc = AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.5e-6
        )

        @test bc.profile == TanhProfile
    end
end

@testset "Type Conversion" begin
    @testset "PeriodicBoundaryCondition Type Conversion" begin
        bc_f32 = PeriodicBoundaryCondition(Float32)
        bc_f64 = FrequencyMaxwell._convert_boundary_type(bc_f32, Float64)

        @test bc_f64 isa PeriodicBoundaryCondition{Float64}
    end

    @testset "AbsorbingBoundaryCondition Type Conversion" begin
        bc_f32 = AbsorbingBoundaryCondition{Float32}(
            Float32(2.0e-6),
            Float32(1.5e-6),
            Float32(0.9),
            TanhProfile
        )

        bc_f64 = FrequencyMaxwell._convert_boundary_type(bc_f32, Float64)

        @test bc_f64 isa AbsorbingBoundaryCondition{Float64}
        @test eltype(bc_f64.thickness) == Float64
        @test eltype(bc_f64.attenuation_thickness) == Float64
        @test eltype(bc_f64.sharpness) == Float64
        @test bc_f64.profile == TanhProfile
    end

    @testset "Mixed Precision Constructor" begin
        # Test that mixed precision inputs are promoted
        bc = AbsorbingBoundaryCondition(
            thickness = 2.0f0,  # Float32
            attenuation_thickness = 1.5,  # Float64
            sharpness = 0.9  # Float64
        )

        # Should promote to Float64
        @test eltype(bc.thickness) == Float64
        @test eltype(bc.attenuation_thickness) == Float64
        @test eltype(bc.sharpness) == Float64
    end
end

@testset "Display Methods" begin
    @testset "PeriodicBoundaryCondition Show" begin
        bc = PeriodicBoundaryCondition()
        output = sprint(show, bc)

        @test occursin("PeriodicBoundaryCondition", output)
        @test occursin("Float64", output)
    end

    @testset "AbsorbingBoundaryCondition Show" begin
        bc = AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.5e-6,
            sharpness = 0.9,
            profile = ExponentialProfile
        )

        output = sprint(show, bc)

        @test occursin("AbsorbingBoundaryCondition", output)
        @test occursin("thickness", output)
        @test occursin("2.0e-6", output)
        @test occursin("sharpness", output)
        @test occursin("0.9", output)
        @test occursin("ExponentialProfile", output)
        @test occursin("subpixel_shift", output)
        @test occursin("0.25", output)
    end
end

@testset "Boundary Condition Integration with Solver" begin
    @testset "Uniform Boundary Conditions" begin
        bc = PeriodicBoundaryCondition()

        solver = ConvergentBornSolver(
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 16),
            boundary_conditions = bc
        )

        # Test that boundary conditions are applied to all dimensions
        @test length(solver.boundary_conditions) == 3
        @test all(bc -> bc isa PeriodicBoundaryCondition, solver.boundary_conditions)
    end

    @testset "Mixed Boundary Conditions" begin
        bc_periodic = PeriodicBoundaryCondition()
        bc_absorbing = AbsorbingBoundaryCondition(
            thickness = 2.0e-6,
            attenuation_thickness = 1.5e-6
        )

        solver = ConvergentBornSolver(
            resolution = (50e-9, 50e-9, 50e-9),
            grid_size = (32, 32, 16),
            boundary_conditions = (bc_periodic, bc_periodic, bc_absorbing)
        )

        # Test per-dimension boundary conditions
        @test solver.boundary_conditions[1] isa PeriodicBoundaryCondition
        @test solver.boundary_conditions[2] isa PeriodicBoundaryCondition
        @test solver.boundary_conditions[3] isa AbsorbingBoundaryCondition
    end
end
