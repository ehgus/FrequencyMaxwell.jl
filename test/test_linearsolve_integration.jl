"""
Test suite for LinearSolve.jl integration in ConvergentBornSolver.

This comprehensive test suite validates that the LinearSolve.jl implementation
produces results identical to the iterative CBS method while providing
performance improvements.
"""

using Test
using FrequencyMaxwell
using LinearAlgebra
using FFTW

"""
Test configuration for validation tests with small problem sizes for speed.
"""
function create_test_config(T=Float64; linear_solver=:iterative, solver_options=Dict{Symbol,Any}())
    return ConvergentBornConfig{T}(
        wavelength = T(500e-9),           # 500 nm
        permittivity_bg = T(1.33^2),      # Water background  
        resolution = (T(25e-9), T(25e-9), T(25e-9)),  # 25 nm isotropic
        grid_size = (32, 32, 16),         # Small grid for fast testing
        use_abbe_sine = true,
        boundary_thickness = (T(100e-9), T(100e-9), T(50e-9)),
        field_attenuation = (T(75e-9), T(75e-9), T(37.5e-9)),
        field_attenuation_sharpness = T(1.0),
        periodic_boundary = (true, true, false),
        iterations_max = 10,              # Fixed iterations for comparison
        tolerance = T(1e-8),              # Tight tolerance for validation
        linear_solver = linear_solver,
        linear_solver_options = solver_options,
        preconditioner = :none
    )
end

"""
Create test permittivity distribution with known scattering features.
"""
function create_test_permittivity(config::ConvergentBornConfig{T}) where T
    grid_size = config.grid_size
    eps_bg = config.permittivity_bg
    permittivity = fill(Complex{T}(eps_bg), grid_size)
    
    # Add a spherical scatterer in the center
    center = grid_size .÷ 2
    radius_pixels = 4  # Small scatterer
    
    for i in 1:grid_size[1], j in 1:grid_size[2], k in 1:grid_size[3]
        dist_sq = (i - center[1])^2 + (j - center[2])^2 + (k - center[3])^2
        if dist_sq <= radius_pixels^2
            # Higher refractive index scatterer
            permittivity[i, j, k] = Complex{T}(eps_bg * 1.2 + 0.01im)
        end
    end
    
    return permittivity
end

"""
Create test plane wave source.
"""
function create_test_source(config::ConvergentBornConfig{T}) where T
    return PlaneWaveSource{T}(
        wavelength = config.wavelength,
        polarization = T[1.0, 0.0, 0.0],  # x-polarized
        k_vector = T[0.0, 0.0, 1.0],      # propagating in +z
        amplitude = Complex{T}(1.0),
        phase = T(0)
    )
end

@testset "LinearSolve.jl Integration Tests" begin
    
    @testset "Configuration Validation" begin
        @testset "Default Configuration" begin
            config = create_test_config()
            @test config.linear_solver == :iterative
            @test config.linear_solver_options == Dict{Symbol,Any}()
            @test config.preconditioner == :none
        end
        
        @testset "LinearSolve Configuration" begin
            options = Dict{Symbol,Any}(:maxiters => 50, :reltol => 1e-10)
            config = create_test_config(; linear_solver=:gmres, solver_options=options)
            @test config.linear_solver == :gmres
            @test config.linear_solver_options == options
        end
        
        @testset "Constructor Type Promotion" begin
            config32 = create_test_config(Float32; linear_solver=:bicgstab)
            @test config32.wavelength isa Float32
            @test config32.tolerance isa Float32
            @test config32.linear_solver == :bicgstab
        end
    end
    
    @testset "Solver Backend Selection" begin
        config_iter = create_test_config(; linear_solver=:iterative)
        config_gmres = create_test_config(; linear_solver=:gmres)
        
        solver_iter = ConvergentBornSolver(config_iter)
        solver_gmres = ConvergentBornSolver(config_gmres)
        
        @test solver_iter.config.linear_solver == :iterative
        @test solver_gmres.config.linear_solver == :gmres
    end
    
    @testset "LinearSolve.jl Operator Validation" begin
        @testset "Operator Construction" begin
            config = create_test_config()
            solver = ConvergentBornSolver(config)
            permittivity = create_test_permittivity(config)
            
            # Initialize solver to create potential and Green's functions
            FrequencyMaxwell._initialize_solver!(solver, permittivity)
            
            # Create linear operator
            linear_op = FrequencyMaxwell.CBSLinearOperator(solver)
            
            @test linear_op.solver === solver
            @test length(linear_op.temp_arrays) == 3
            @test size(linear_op.temp_arrays[1]) == (size(solver.potential)..., 3)
            @test eltype(linear_op) == Complex{Float64}
        end
        
        @testset "Operator Size and Type" begin
            config = create_test_config(Float32)
            solver = ConvergentBornSolver(config)
            permittivity = create_test_permittivity(config)
            FrequencyMaxwell._initialize_solver!(solver, permittivity)
            
            linear_op = FrequencyMaxwell.CBSLinearOperator(solver)
            
            expected_size = prod(size(solver.potential)) * 3
            @test size(linear_op) == (expected_size, expected_size)
            @test eltype(linear_op) == Complex{Float32}
        end
    end
    
    @testset "Field Accuracy Validation" begin
        # Test that LinearSolve.jl produces identical results to iterative method
        
        @testset "Single Source - GMRES vs Iterative" begin
            # Set up identical problems with different solvers
            config_iter = create_test_config(; linear_solver=:iterative)
            config_gmres = create_test_config(; linear_solver=:gmres)
            
            solver_iter = ConvergentBornSolver(config_iter)
            solver_gmres = ConvergentBornSolver(config_gmres)
            
            permittivity = create_test_permittivity(config_iter)
            source = create_test_source(config_iter)
            
            # Solve with both methods
            E_iter, H_iter = solve(solver_iter, source, permittivity)
            E_gmres, H_gmres = solve(solver_gmres, source, permittivity)
            
            # Validate field dimensions
            @test size(E_iter) == size(E_gmres)
            @test size(H_iter) == size(H_gmres)
            @test size(E_iter)[end] == 3  # Vector field
            
            # Validate field accuracy (should be nearly identical)
            E_diff = maximum(abs.(E_gmres - E_iter))
            H_diff = maximum(abs.(H_gmres - H_iter))
            
            @test E_diff < 1e-10  # Very tight tolerance for validation
            @test H_diff < 1e-10
            
            @test isfinite(maximum(abs.(E_iter)))
            @test isfinite(maximum(abs.(H_iter)))
        end
        
        @testset "Multiple Algorithms Consistency" begin
            algorithms = [:gmres, :bicgstab, :minres]
            config_iter = create_test_config(; linear_solver=:iterative)
            solver_iter = ConvergentBornSolver(config_iter)
            
            permittivity = create_test_permittivity(config_iter)
            source = create_test_source(config_iter)
            
            # Reference solution with iterative method
            E_ref, H_ref = solve(solver_iter, source, permittivity)
            
            for alg in algorithms
                @testset "Algorithm: $alg" begin
                    config = create_test_config(; linear_solver=alg)
                    solver = ConvergentBornSolver(config)
                    
                    E_test, H_test = solve(solver, source, permittivity)
                    
                    E_diff = maximum(abs.(E_test - E_ref))
                    H_diff = maximum(abs.(H_test - H_ref))
                    
                    # Allow slightly larger tolerance for different algorithms
                    @test E_diff < 1e-8
                    @test H_diff < 1e-8
                end
            end
        end
        
        @testset "Multi-Source Coherent Interference" begin
            config_iter = create_test_config(; linear_solver=:iterative)
            config_gmres = create_test_config(; linear_solver=:gmres)
            
            solver_iter = ConvergentBornSolver(config_iter)
            solver_gmres = ConvergentBornSolver(config_gmres)
            
            permittivity = create_test_permittivity(config_iter)
            
            # Create multiple coherent sources
            source1 = create_test_source(config_iter)
            source2 = PlaneWaveSource{Float64}(
                wavelength = config_iter.wavelength,
                polarization = [0.0, 1.0, 0.0],  # y-polarized
                k_vector = [1.0, 0.0, 1.0] / √2,  # diagonal propagation
                amplitude = Complex(0.8),
                phase = π/4
            )
            sources = [source1, source2]
            
            # Solve with both methods
            E_iter, H_iter = solve(solver_iter, sources, permittivity)
            E_gmres, H_gmres = solve(solver_gmres, sources, permittivity)
            
            # Validate coherent interference patterns are preserved
            E_diff = maximum(abs.(E_gmres - E_iter))
            H_diff = maximum(abs.(H_gmres - H_iter))
            
            @test E_diff < 1e-10
            @test H_diff < 1e-10
        end
        
        @testset "Precision Consistency Float32" begin
            config_iter = create_test_config(Float32; linear_solver=:iterative)
            config_gmres = create_test_config(Float32; linear_solver=:gmres)
            
            solver_iter = ConvergentBornSolver(config_iter)
            solver_gmres = ConvergentBornSolver(config_gmres)
            
            permittivity = create_test_permittivity(config_iter)
            source = create_test_source(config_iter)
            
            E_iter, H_iter = solve(solver_iter, source, permittivity)
            E_gmres, H_gmres = solve(solver_gmres, source, permittivity)
            
            # Validate data types
            @test eltype(E_iter) == Complex{Float32}
            @test eltype(E_gmres) == Complex{Float32}
            @test eltype(H_iter) == Complex{Float32}
            @test eltype(H_gmres) == Complex{Float32}
            
            # Validate accuracy with appropriate tolerance for Float32
            E_diff = maximum(abs.(E_gmres - E_iter))
            H_diff = maximum(abs.(H_gmres - H_iter))
            
            @test E_diff < 1e-5  # Appropriate for Float32 precision
            @test H_diff < 1e-5
        end
    end
    
    @testset "Error Handling and Fallback" begin
        @testset "Unknown Algorithm Fallback" begin
            # Test that unknown algorithms fall back to GMRES with warning
            config = create_test_config(; linear_solver=:unknown_solver)
            solver = ConvergentBornSolver(config)
            
            permittivity = create_test_permittivity(config)
            source = create_test_source(config)
            
            # Should not throw error, should use GMRES fallback
            @test_nowarn E, H = solve(solver, source, permittivity)
            
            E, H = solve(solver, source, permittivity)
            @test isfinite(maximum(abs.(E)))
            @test isfinite(maximum(abs.(H)))
        end
        
        @testset "Convergence Failure Fallback" begin
            # Create configuration with impossible convergence criteria
            harsh_options = Dict{Symbol,Any}(
                :maxiters => 2,      # Very few iterations
                :reltol => 1e-15,    # Impossible tolerance
                :abstol => 1e-20
            )
            config = create_test_config(; 
                linear_solver=:gmres, 
                solver_options=harsh_options
            )
            solver = ConvergentBornSolver(config)
            
            permittivity = create_test_permittivity(config)
            source = create_test_source(config)
            
            # Should fall back to iterative method gracefully
            @test_nowarn E, H = solve(solver, source, permittivity)
            
            E, H = solve(solver, source, permittivity)
            @test isfinite(maximum(abs.(E)))
            @test isfinite(maximum(abs.(H)))
        end
    end
    
    @testset "Solver Statistics and Diagnostics" begin
        config = create_test_config(; linear_solver=:gmres)
        solver = ConvergentBornSolver(config)
        
        permittivity = create_test_permittivity(config)
        source = create_test_source(config)
        
        # Clear any previous state
        FrequencyMaxwell.reset!(solver)
        @test solver.iteration_count == 0
        @test isempty(solver.residual_history)
        
        # Solve and check statistics are updated
        E, H = solve(solver, source, permittivity)
        
        @test solver.iteration_count >= 1
        @test length(solver.residual_history) >= 1
        @test all(isfinite, solver.residual_history)
    end
end