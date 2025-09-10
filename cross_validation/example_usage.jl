#!/usr/bin/env julia

"""
Cross-Validation Framework Example Usage

This example demonstrates how to use the FrequencyMaxwell cross-validation 
framework with automatic GitHub repository integration.

The framework automatically downloads and uses the Helmholtz adjoint solver
from: https://github.com/ehgus/Helmholtz-adjoint-solver

Prerequisites:
1. Julia (>= 1.8)
2. MATLAB with MATLAB.jl package
3. Git (for repository cloning)

Setup:
    julia setup.jl --clone-repo --install-deps

Usage:
    julia example_usage.jl
"""

# Add current directory to path for local modules
push!(LOAD_PATH, ".")

using Pkg

# Import required packages 
# If packages are missing, run: julia -e "using Pkg; Pkg.add([\"MATLAB\", \"JSON3\", \"Dates\"])"
using MATLAB, JSON3, Dates, LinearAlgebra, Statistics

# Load the cross-validation framework
include("src/CrossValidation.jl")
using .CrossValidation

println("="^70)
println("FrequencyMaxwell Cross-Validation Framework")
println("GitHub Repository: https://github.com/ehgus/Helmholtz-adjoint-solver")
println("="^70)

"""
    run_basic_validation_example()

Run a basic validation example with a few test cases.
"""
function run_basic_validation_example()
    println("\n🚀 Setting up cross-validation framework...")

    # Initialize framework (automatically uses GitHub repository)
    cv = CrossValidationFramework(
        reference_data_path = "./reference_data",
        reports_path = "./reports"
    )

    println("✅ Framework initialized with GitHub repository integration")

    # Add forward validation test cases
    println("\n📋 Adding forward validation test cases...")

    add_test_case!(cv,
        TestCase(
            "metalens_forward",
            "forward",
            Dict{String, Any}(
                "wavelength" => 532e-9,
                "focal_length" => 50e-6,
                "grid_size" => [32.0, 32.0]  # Small for demo
            ),
            Dict{String, Float64}(
                "field_relative_error" => 1e-4,
                "energy_relative_error" => 1e-3,
                "convergence_threshold" => 1e-8
            ),
            "Basic metalens forward solver validation with 32x32 grid"
        ))

    add_test_case!(cv,
        TestCase(
            "grating_forward",
            "forward",
            Dict{String, Any}(
                "wavelength" => 532e-9,
                "grating_period" => 2e-6,
                "grid_size" => [64.0, 16.0]  # Optimized for grating
            ),
            Dict{String, Float64}(
                "field_relative_error" => 1e-4,
                "energy_relative_error" => 1e-3
            ),
            "Basic diffraction grating forward solver validation"
        ))

    # Add adjoint validation test case
    println("📋 Adding adjoint validation test cases...")

    add_test_case!(cv,
        TestCase(
            "metalens_adjoint",
            "adjoint",
            Dict{String, Any}(
                "wavelength" => 532e-9,
                "target_focal_length" => 50e-6,
                "grid_size" => [24.0, 24.0],  # Smaller for adjoint
                "max_iterations" => 5     # Quick demo
            ),
            Dict{String, Float64}(
                "gradient_relative_error" => 1e-3,
                "convergence_threshold" => 1e-6
            ),
            "Basic metalens adjoint optimization validation"
        ))

    println("✅ Added $(length(cv.test_cases)) test cases")

    # Run validation
    println("\n🧪 Running cross-validation tests...")
    println("This will:")
    println("  1. Initialize MATLAB session")
    println("  2. Load solver from GitHub repository")
    println("  3. Execute test cases")
    println("  4. Generate validation reports")

    try
        results = run_validation!(cv,
            generate_reports = true,
            save_reference = false  # Don't overwrite reference data
        )

        # Print summary
        println("\n" * "="^50)
        println("VALIDATION SUMMARY")
        println("="^50)

        passed = sum(r.success for r in results)
        total = length(results)
        success_rate = round(100 * passed / total, digits = 1)

        println("Tests run: $total")
        println("Passed: $passed")
        println("Failed: $(total - passed)")
        println("Success rate: $success_rate%")

        if passed == total
            println("\n🎉 All tests passed!")
        else
            println("\n⚠️  Some tests failed. Check the detailed report.")
        end

        # Show report location
        report_files = filter(f -> endswith(f, ".md"), readdir("./reports", join = true))
        if !isempty(report_files)
            latest_report = sort(report_files, by = mtime, rev = true)[1]
            println("\n📄 Detailed report: $latest_report")
        end

        return results

    catch e
        println("\n❌ Validation failed with error:")
        println("Error: $e")

        # Check common issues
        println("\n🔍 Troubleshooting:")

        if !isdir("./Helmholtz-adjoint-solver")
            println("  • Repository not found. Run: julia setup.jl --clone-repo")
        end

        # MATLAB.jl is available (successfully imported at top level)
        println("  • MATLAB.jl is available")

        rethrow(e)
    end
end

"""
    run_comprehensive_validation_example()

Run a comprehensive validation with all available test cases.
"""
function run_comprehensive_validation_example()
    println("\n🚀 Setting up comprehensive validation...")

    cv = CrossValidationFramework()

    # Add all available test cases
    test_cases = [
        # Forward tests
        ("metalens_forward", "forward",
            Dict(
                "wavelength" => 532e-9,
                "focal_length" => 50e-6,
                "grid_size" => [64, 64]
            ),
            "Metalens forward solver validation"), ("grating_forward",
            "forward",
            Dict(
                "wavelength" => 632e-9,
                "grating_period" => 1.5e-6,
                "grid_size" => [128, 32]
            ),
            "Grating forward solver validation"),
        ("two_beam_forward", "forward",
            Dict(
                "wavelength" => 532e-9,
                "beam_angle" => 15.0,
                "grid_size" => [64, 64]
            ),
            "Two-beam interference validation"),
        ("sio2_sphere_forward",
            "forward",
            Dict(
                "wavelength" => 532e-9,
                "sphere_radius" => 2.5e-6,
                "grid_size" => [64, 64]
            ),
            "SiO2 sphere scattering validation"),
        ("helical_metalens_forward",
            "forward",
            Dict(
                "wavelength" => 532e-9,
                "topological_charge" => 1,
                "focal_length" => 50e-6,
                "grid_size" => [64, 64]
            ),
            "Helical metalens forward validation"),

        # Adjoint tests
        ("metalens_adjoint",
            "adjoint",
            Dict(
                "wavelength" => 532e-9,
                "target_focal_length" => 50e-6,
                "grid_size" => [32, 32],
                "max_iterations" => 10
            ),
            "Metalens adjoint optimization"), ("grating_adjoint",
            "adjoint",
            Dict(
                "wavelength" => 632e-9,
                "target_efficiency" => 0.8,
                "grid_size" => [64, 16],
                "max_iterations" => 15
            ),
            "Grating adjoint optimization"),
        ("double_helix_adjoint",
            "adjoint",
            Dict(
                "wavelength" => 532e-9,
                "helix_separation" => 5e-6,
                "grid_size" => [48, 48],
                "max_iterations" => 8
            ),
            "Double helix PSF optimization"),
        ("single_helix_adjoint",
            "adjoint",
            Dict(
                "wavelength" => 532e-9,
                "topological_charge" => 2,
                "grid_size" => [40, 40],
                "max_iterations" => 12
            ),
            "Single helix PSF optimization")
    ]

    for (name, test_type, params, description) in test_cases
        # Set appropriate tolerances based on test type
        if test_type == "forward"
            tolerances = Dict{String, Float64}(
                "field_relative_error" => 1e-6,
                "energy_relative_error" => 1e-8,
                "convergence_threshold" => 1e-10
            )
        else  # adjoint
            tolerances = Dict{String, Float64}(
                "gradient_relative_error" => 1e-5,
                "convergence_threshold" => 1e-8
            )
        end

        add_test_case!(cv, TestCase(name, test_type, params, tolerances, description))
    end

    println("✅ Added $(length(cv.test_cases)) comprehensive test cases")

    # Run validation
    results = run_validation!(cv,
        generate_reports = true,
        save_reference = true  # Save as new reference
    )

    return results
end

"""
    demonstrate_repository_management()

Demonstrate GitHub repository management features.
"""
function demonstrate_repository_management()
    println("\n📂 Repository Management Demo")
    println("="^40)

    # Check if repository exists
    repo_path = "./Helmholtz-adjoint-solver"

    if isdir(repo_path)
        println("✅ Repository found at: $repo_path")

        # Show repository information
        cd(repo_path) do
            remote_url = strip(read(`git remote get-url origin`, String))
            current_branch = strip(read(`git branch --show-current`, String))
            last_commit = strip(read(`git log -1 --format="%h %s"`, String))

            println("📍 Remote: $remote_url")
            println("🌿 Branch: $current_branch")
            println("💾 Last commit: $last_commit")
        end

        # List key files
        println("\n📁 Repository contents:")
        for item in readdir(repo_path)[1:min(10, length(readdir(repo_path)))]
            item_path = joinpath(repo_path, item)
            icon = isdir(item_path) ? "📁" : "📄"
            println("  $icon $item")
        end

        if length(readdir(repo_path)) > 10
            println("  ... and $(length(readdir(repo_path)) - 10) more items")
        end

    else
        println("❌ Repository not found at: $repo_path")
        println("Run: julia setup.jl --clone-repo")
    end
end

"""
    main()

Main function to run the example.
"""
function main()
    println("FrequencyMaxwell Cross-Validation Framework Example")
    println("Using GitHub repository: https://github.com/ehgus/Helmholtz-adjoint-solver")

    # Show repository status
    demonstrate_repository_management()

    # Ask user what to run
    println("\nSelect validation type:")
    println("1. Basic validation (3 test cases, ~2-5 minutes)")
    println("2. Comprehensive validation (all test cases, ~10-20 minutes)")
    println("3. Repository management demo only")
    print("Choice (1/2/3): ")

    choice = readline()

    if choice == "1"
        println("\n🎯 Running basic validation example...")
        results = run_basic_validation_example()

    elseif choice == "2"
        println("\n🎯 Running comprehensive validation example...")
        results = run_comprehensive_validation_example()

    elseif choice == "3"
        println("\n📂 Repository management demo completed.")
        return

    else
        println("Invalid choice. Running basic validation by default...")
        results = run_basic_validation_example()
    end

    # Final summary
    if @isdefined(results)
        println("\n" * "="^70)
        println("FINAL SUMMARY")
        println("="^70)
        println("Framework successfully used GitHub repository:")
        println("  https://github.com/ehgus/Helmholtz-adjoint-solver")
        println("")
        println("Validation results:")
        println("  Total tests: $(length(results))")
        println("  Passed: $(sum(r.success for r in results))")
        println("  Success rate: $(round(100 * sum(r.success for r in results) / length(results), digits=1))%")

        # Show next steps
        println("\n🔄 Next steps:")
        println("  • Check detailed reports in ./reports/")
        println("  • Review reference data in ./reference_data/")
        println("  • Update repository: cd Helmholtz-adjoint-solver && git pull")
        println("  • Customize test cases by editing src/test_cases.jl")
    end
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
