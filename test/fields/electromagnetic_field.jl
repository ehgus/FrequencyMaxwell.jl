"""
ElectromagneticField structure tests.
Tests field construction, utilities, and domain properties.
"""

@testset "ElectromagneticField Structure" begin
    @testset "Field Construction" begin
        # Test basic field construction with correct constructor
        Efield = rand(ComplexF64, 32, 32, 16, 3)  # Grid + 3 components
        Hfield = rand(ComplexF64, 32, 32, 16, 3)  # Grid + 3 components
        grid_size = (32, 32, 16)
        resolution = (100e-9, 100e-9, 200e-9)
        wavelength = 500e-9

        EMfield = ElectromagneticField(Efield, Hfield, grid_size, resolution, wavelength)

        @test size(EMfield.E) == (32, 32, 16, 3)
        @test size(EMfield.H) == (32, 32, 16, 3)
        @test EMfield.grid_size == grid_size
        @test EMfield.resolution == resolution
        @test EMfield.wavelength ≈ wavelength
    end

    @testset "Field Utilities" begin
        # Create test fields with correct constructor
        Efield = ones(ComplexF64, 16, 16, 8, 3) * 0.1
        Hfield = ones(ComplexF64, 16, 16, 8, 3) * 0.01
        grid_size = (16, 16, 8)
        resolution = (50e-9, 50e-9, 100e-9)
        wavelength = 633e-9

        EMfield = ElectromagneticField(Efield, Hfield, grid_size, resolution, wavelength)

        # Test energy calculation
        energy = field_energy(EMfield)
        @test energy > 0
        @test isfinite(energy)

        # Test intensity calculation
        intensity = field_intensity(EMfield)
        @test size(intensity) == (16, 16, 8)
        @test all(intensity .≥ 0)

        # Test Poynting vector
        S = poynting_vector(EMfield)
        @test size(S) == (16, 16, 8, 3)
    end

    @testset "Field Domain Properties" begin
        Efield = rand(ComplexF64, 20, 30, 10, 3)
        Hfield = rand(ComplexF64, 20, 30, 10, 3)
        grid_size = (20, 30, 10)
        resolution = (25e-9, 40e-9, 80e-9)
        wavelength = 500e-9

        EMfield = ElectromagneticField(Efield, Hfield, grid_size, resolution, wavelength)

        domain = domain_size(EMfield)
        @test domain[1] ≈ 20 * 25e-9
        @test domain[2] ≈ 30 * 40e-9
        @test domain[3] ≈ 10 * 80e-9
    end
end
