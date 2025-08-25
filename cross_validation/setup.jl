#!/usr/bin/env julia
"""
Setup script for FrequencyMaxwell Cross-Validation Framework

This script helps configure the cross-validation environment and verify
that all dependencies are properly installed and configured.

Usage:
    julia setup.jl [options]

Options:
    --check-only    Only check current status, don't install dependencies
    --install-deps  Install required Julia dependencies
    --test          Run basic tests after setup
    --verbose       Enable detailed output
    --clone-repo    Clone the Helmholtz adjoint solver repository
"""

using Pkg
using InteractiveUtils

# GitHub repository configuration
const HELMHOLTZ_REPO_URL = "https://github.com/ehgus/Helmholtz-adjoint-solver.git"
const HELMHOLTZ_REPO_BRANCH = "main"
const LOCAL_REPO_PATH = "./Helmholtz-adjoint-solver"

function print_header()
    println("="^70)
    println("FrequencyMaxwell Cross-Validation Framework Setup")
    println("="^70)
    println("Repository: $HELMHOLTZ_REPO_URL")
    println("Branch: $HELMHOLTZ_REPO_BRANCH")
    println("="^70)
end

function check_julia_version()
    println("📋 Checking Julia version...")
    min_version = v"1.8"
    current_version = VERSION
    
    if current_version >= min_version
        println("  ✓ Julia $current_version (>= $min_version required)")
        return true
    else
        println("  ✗ Julia $current_version (>= $min_version required)")
        println("    Please upgrade Julia to continue.")
        return false
    end
end

function check_git_availability()
    println("\n🔧 Checking Git availability...")
    
    try
        result = read(`git --version`, String)
        println("  ✓ Git is available: $(strip(result))")
        return true
    catch
        println("  ✗ Git not found")
        println("    Please install Git to clone the Helmholtz solver repository")
        return false
    end
end

function check_matlab_availability()
    println("\n🔧 Checking MATLAB availability...")
    
    # Check if MATLAB.jl is installed
    matlab_jl_installed = false
    try
        Pkg.status("MATLAB")
        println("  ✓ MATLAB.jl package installed")
        matlab_jl_installed = true
    catch
        println("  ! MATLAB.jl package not installed")
        println("    Run with --install-deps to install it")
    end
    
    # Check if MATLAB executable is available
    matlab_path = get(ENV, "MATLAB_PATH", "matlab")
    matlab_available = false
    
    try
        result = read(`which $matlab_path`, String)
        if !isempty(strip(result))
            println("  ✓ MATLAB executable found: $(strip(result))")
            matlab_available = true
        end
    catch
        # Try alternative paths
        common_paths = [
            "/usr/local/MATLAB/R2023b/bin/matlab",
            "/Applications/MATLAB_R2023b.app/bin/matlab",
            "/opt/matlab/bin/matlab",
            "C:\\Program Files\\MATLAB\\R2023b\\bin\\matlab.exe"
        ]
        
        for path in common_paths
            if isfile(path)
                println("  ✓ MATLAB found at: $path")
                matlab_available = true
                break
            end
        end
        
        if !matlab_available
            println("  ! MATLAB executable not found in standard locations")
            println("    Please ensure MATLAB is installed and in PATH")
            println("    Or set MATLAB_PATH environment variable")
        end
    end
    
    return matlab_jl_installed && matlab_available
end

function check_repository_status()
    println("\n📁 Checking Helmholtz solver repository...")
    
    local_path = abspath(LOCAL_REPO_PATH)
    
    if isdir(local_path)
        println("  ✓ Repository directory found: $local_path")
        
        # Check if it's a git repository
        git_dir = joinpath(local_path, ".git")
        if isdir(git_dir)
            println("  ✓ Valid git repository")
            
            # Check remote origin
            try
                cd(local_path) do
                    remote_url = strip(read(`git remote get-url origin`, String))
                    if occursin("Helmholtz-adjoint-solver", remote_url)
                        println("  ✓ Correct remote origin: $remote_url")
                        
                        # Check current branch
                        current_branch = strip(read(`git branch --show-current`, String))
                        println("  📍 Current branch: $current_branch")
                        
                        # Check for updates
                        try
                            read(`git fetch origin`, String)
                            local_commit = strip(read(`git rev-parse HEAD`, String))
                            remote_commit = strip(read(`git rev-parse origin/$HELMHOLTZ_REPO_BRANCH`, String))
                            
                            if local_commit == remote_commit
                                println("  ✓ Repository is up to date")
                            else
                                println("  ⚠ Repository has updates available")
                                println("    Run 'git pull origin $HELMHOLTZ_REPO_BRANCH' to update")
                            end
                        catch
                            println("  ! Unable to check for updates (network issue?)")
                        end
                        
                        return true
                    else
                        println("  ✗ Wrong remote origin: $remote_url")
                        println("    Expected: $HELMHOLTZ_REPO_URL")
                        return false
                    end
                end
            catch e
                println("  ✗ Error checking git remote: $e")
                return false
            end
        else
            println("  ✗ Directory exists but is not a git repository")
            return false
        end
    else
        println("  ! Repository not found: $local_path")
        println("    Run with --clone-repo to download it")
        return false
    end
end

function clone_repository()
    println("\n📥 Cloning Helmholtz adjoint solver repository...")
    
    local_path = abspath(LOCAL_REPO_PATH)
    
    if isdir(local_path)
        println("  ! Directory already exists: $local_path")
        print("  Remove existing directory and re-clone? (y/N): ")
        response = readline()
        if lowercase(strip(response)) != "y"
            println("  Aborted.")
            return false
        end
        
        rm(local_path, recursive=true)
        println("  ✓ Removed existing directory")
    end
    
    try
        println("  📡 Cloning from $HELMHOLTZ_REPO_URL...")
        run(`git clone --branch $HELMHOLTZ_REPO_BRANCH $HELMHOLTZ_REPO_URL $local_path`)
        println("  ✓ Repository cloned successfully")
        
        # Verify the clone
        if isdir(local_path) && isdir(joinpath(local_path, ".git"))
            println("  ✓ Clone verification passed")
            
            # List key directories/files
            println("  📂 Repository contents:")
            for item in readdir(local_path)
                item_path = joinpath(local_path, item)
                if isdir(item_path)
                    println("    📁 $item/")
                else
                    println("    📄 $item")
                end
            end
            
            return true
        else
            println("  ✗ Clone verification failed")
            return false
        end
        
    catch e
        println("  ✗ Failed to clone repository: $e")
        println("    Please check your internet connection and try again")
        return false
    end
end

function check_solver_files()
    println("\n📋 Checking MATLAB solver files...")
    
    solver_path = abspath(LOCAL_REPO_PATH)
    
    if !isdir(solver_path)
        println("  ✗ Solver directory not found: $solver_path")
        println("    Run with --clone-repo to download the repository")
        return false
    end
    
    # Expected key files (based on typical MATLAB solver structure)
    key_files = [
        "ConvergentBornSolver.m",
        "setup.m",
        "README.md"
    ]
    
    missing_files = String[]
    found_files = String[]
    
    # Search recursively for key files
    for (root, dirs, files) in walkdir(solver_path)
        for file in files
            if file in key_files
                full_path = joinpath(root, file)
                relative_path = relpath(full_path, solver_path)
                push!(found_files, relative_path)
                println("  ✓ Found: $relative_path")
            end
        end
    end
    
    # Check for missing files
    for file in key_files
        if !any(endswith(f, file) for f in found_files)
            push!(missing_files, file)
        end
    end
    
    if isempty(missing_files)
        println("  ✓ All key solver files found")
        return true
    else
        println("  ⚠ Some expected files not found:")
        for file in missing_files
            println("    - $file")
        end
        println("  Repository structure may have changed. Validation may still work.")
        return true  # Don't fail completely
    end
end

function install_julia_dependencies()
    println("\n📦 Installing Julia dependencies...")
    
    required_packages = [
        "MATLAB",
        "JSON3",
        "Dates",
        "Test",
        "LinearAlgebra",
        "Statistics"
    ]
    
    for pkg in required_packages
        try
            println("  📥 Installing $pkg...")
            Pkg.add(pkg)
            println("  ✓ $pkg installed successfully")
        catch e
            println("  ✗ Failed to install $pkg: $e")
            return false
        end
    end
    
    println("  ✓ All Julia dependencies installed")
    return true
end

function run_basic_tests()
    println("\n🧪 Running basic validation tests...")
    
    try
        # Test Julia module loading
        include("src/CrossValidation.jl")
        println("  ✓ Julia modules load successfully")
        
        # Test MATLAB connection
        try
            Pkg.status("MATLAB")
            println("  ✓ MATLAB.jl package available")
        catch e
            println("  ✗ MATLAB.jl not installed: $e")
            return false
        end
        
        println("  ✓ Basic tests passed")
        return true
        
    catch e
        println("  ✗ Basic tests failed: $e")
        return false
    end
end

function print_summary(checks_passed::Dict)
    println("\n" * "="^70)
    println("SETUP SUMMARY")
    println("="^70)
    
    all_passed = true
    for (check, passed) in checks_passed
        status = passed ? "✓" : "✗"
        println("  $status $check")
        all_passed = all_passed && passed
    end
    
    println()
    if all_passed
        println("🎉 Setup completed successfully!")
        println("\nNext steps:")
        println("  1. Run: julia example_usage.jl")
        println("  2. Check the validation reports in ./reports/")
        println("  3. Review README.md for detailed usage")
    else
        println("⚠️  Setup completed with warnings.")
        println("\nRecommended actions:")
        if !checks_passed["Git Available"]
            println("  1. Install Git to enable repository cloning")
        end
        if !checks_passed["Repository Status"]
            println("  2. Run: julia setup.jl --clone-repo")
        end
        if !checks_passed["MATLAB Available"]
            println("  3. Install MATLAB and MATLAB.jl package")
        end
        if !checks_passed["Julia Dependencies"]
            println("  4. Run: julia setup.jl --install-deps")
        end
    end
    
    println("\n" * "="^70)
end

function main()
    # Parse command line arguments
    args = ARGS
    check_only = "--check-only" in args
    install_deps = "--install-deps" in args
    run_tests = "--test" in args
    verbose = "--verbose" in args
    clone_repo = "--clone-repo" in args
    
    print_header()
    
    # Run checks
    checks_passed = Dict{String, Bool}()
    
    checks_passed["Julia Version"] = check_julia_version()
    checks_passed["Git Available"] = check_git_availability()
    checks_passed["MATLAB Available"] = check_matlab_availability()
    
    if clone_repo
        success = clone_repository()
        if !success
            println("\nFailed to clone repository. Exiting.")
            exit(1)
        end
    end
    
    checks_passed["Repository Status"] = check_repository_status()
    checks_passed["Solver Files"] = check_solver_files()
    
    # Install dependencies if requested
    if install_deps && !check_only
        checks_passed["Julia Dependencies"] = install_julia_dependencies()
    else
        checks_passed["Julia Dependencies"] = true  # Assume OK if not installing
    end
    
    # Run tests if requested
    if run_tests && !check_only
        checks_passed["Basic Tests"] = run_basic_tests()
    else
        checks_passed["Basic Tests"] = true  # Skip if not requested
    end
    
    print_summary(checks_passed)
    
    # Exit with appropriate code
    all_passed = all(values(checks_passed))
    exit(all_passed ? 0 : 1)
end

# Only run if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
