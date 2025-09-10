"""
Repository Manager Module for Cross-Validation Framework

This module provides robust GitHub repository integration with comprehensive
fallback mechanisms, automatic updates, version management, and offline capabilities.

Features:
- Automatic repository cloning and management
- Version pinning and branch management
- Fallback mechanisms for network issues
- Offline mode with cached repositories
- Integrity verification and validation
- Automatic updates with conflict resolution
- Multiple repository mirror support
"""

module RepositoryManager

using Dates
using Logging
using JSON3

export RepositoryConfig, RepositoryManager, clone_repository, update_repository
export check_repository_integrity, setup_fallback_mirrors, enable_offline_mode
export get_repository_info, sync_repository, validate_repository_structure

"""
Repository configuration structure
"""
struct RepositoryConfig
    primary_url::String
    fallback_urls::Vector{String}
    local_path::String
    branch::String
    commit_hash::Union{String, Nothing}
    auto_update::Bool
    verify_integrity::Bool
    offline_mode::Bool
    cache_path::String

    function RepositoryConfig(;
            primary_url = "https://github.com/ehgus/Helmholtz-adjoint-solver.git",
            fallback_urls = String[],
            local_path = "./Helmholtz-adjoint-solver",
            branch = "main",
            commit_hash = nothing,
            auto_update = true,
            verify_integrity = true,
            offline_mode = false,
            cache_path = "./repository_cache"
    )
        new(primary_url, fallback_urls, local_path, branch, commit_hash,
            auto_update, verify_integrity, offline_mode, cache_path)
    end
end

"""
Repository manager with comprehensive functionality
"""
mutable struct RepositoryManager
    config::RepositoryConfig
    status::Dict{String, Any}
    last_update::Union{DateTime, Nothing}
    integrity_verified::Bool
    fallback_active::Bool
    offline_available::Bool

    function RepositoryManager(config::RepositoryConfig)
        status = Dict{String, Any}(
            "initialized" => false,
            "cloned" => false,
            "up_to_date" => false,
            "integrity_ok" => false
        )

        new(config, status, nothing, false, false, false)
    end
end

"""
    clone_repository(manager::RepositoryManager; force=false)

Clone repository with comprehensive fallback mechanisms.
"""
function clone_repository(manager::RepositoryManager; force = false)
    @info "Cloning repository with fallback support..."

    # Check if repository already exists
    if isdir(manager.config.local_path) && !force
        @info "Repository already exists at $(manager.config.local_path)"
        return update_repository(manager)
    end

    # Remove existing directory if force is true
    if force && isdir(manager.config.local_path)
        @info "Removing existing repository (force=true)"
        rm(manager.config.local_path, recursive = true, force = true)
    end

    # Try primary URL first
    success = try_clone_from_url(manager, manager.config.primary_url)

    if !success && !isempty(manager.config.fallback_urls)
        @warn "Primary repository clone failed, trying fallback URLs..."

        for fallback_url in manager.config.fallback_urls
            @info "Trying fallback URL: $fallback_url"
            success = try_clone_from_url(manager, fallback_url)
            if success
                manager.fallback_active = true
                break
            end
        end
    end

    if !success
        # Try offline mode if available
        if manager.config.offline_mode && check_offline_repository(manager)
            @info "Using offline cached repository"
            success = setup_offline_repository(manager)
        end
    end

    if success
        manager.status["cloned"] = true
        manager.status["initialized"] = true

        # Verify integrity if requested
        if manager.config.verify_integrity
            verify_repository_integrity(manager)
        end

        # Pin to specific commit if specified
        if manager.config.commit_hash !== nothing
            pin_to_commit(manager, manager.config.commit_hash)
        end

        @info "Repository cloned successfully"
        return true
    else
        @error "Failed to clone repository from any source"
        return false
    end
end

"""
    try_clone_from_url(manager::RepositoryManager, url::String)

Attempt to clone repository from a specific URL.
"""
function try_clone_from_url(manager::RepositoryManager, url::String)
    try
        @info "Cloning from: $url"

        # Construct git clone command
        clone_cmd = if manager.config.branch != "main"
            `git clone --branch $(manager.config.branch) $url $(manager.config.local_path)`
        else
            `git clone $url $(manager.config.local_path)`
        end

        # Set timeout for clone operation
        clone_process = run(clone_cmd, wait = false)

        # Wait for completion with timeout (5 minutes)
        timeout = 300  # seconds
        start_time = time()

        while process_running(clone_process) && (time() - start_time) < timeout
            sleep(1)
        end

        if process_running(clone_process)
            kill(clone_process)
            @warn "Clone operation timed out"
            return false
        end

        if clone_process.exitcode == 0
            @info "Clone successful from $url"
            return true
        else
            @warn "Clone failed from $url (exit code: $(clone_process.exitcode))"
            return false
        end

    catch e
        @warn "Error cloning from $url: $e"
        return false
    end
end

"""
    update_repository(manager::RepositoryManager)

Update repository with conflict resolution and fallback support.
"""
function update_repository(manager::RepositoryManager)
    if !isdir(manager.config.local_path)
        @warn "Repository not found, attempting to clone..."
        return clone_repository(manager)
    end

    @info "Updating repository..."

    try
        cd(manager.config.local_path) do
            # Check repository status
            status_result = read(`git status --porcelain`, String)

            if !isempty(strip(status_result))
                @warn "Repository has uncommitted changes:"
                println(status_result)

                # Stash changes
                @info "Stashing local changes..."
                run(`git stash push -m "Auto-stash before update $(now())"`)
            end

            # Fetch latest changes
            @info "Fetching latest changes..."

            # Try primary remote first
            fetch_success = try
                run(`git fetch origin`)
                true
            catch e
                @warn "Failed to fetch from primary remote: $e"
                false
            end

            # Try fallback remotes if primary fails
            if !fetch_success && !isempty(manager.config.fallback_urls)
                fetch_success = try_fallback_fetch(manager)
            end

            if !fetch_success
                @error "Failed to fetch from any remote"
                return false
            end

            # Check for updates
            local_commit = strip(read(`git rev-parse HEAD`, String))
            remote_commit = strip(read(`git rev-parse origin/$(manager.config.branch)`, String))

            if local_commit == remote_commit
                @info "Repository is already up to date"
                manager.status["up_to_date"] = true
                return true
            end

            # Perform update
            @info "Updating to latest version..."

            if manager.config.commit_hash !== nothing
                # Update to specific commit
                run(`git checkout $(manager.config.commit_hash)`)
            else
                # Update to latest on branch
                run(`git merge origin/$(manager.config.branch)`)
            end

            manager.last_update = now()
            manager.status["up_to_date"] = true

            # Verify integrity after update
            if manager.config.verify_integrity
                verify_repository_integrity(manager)
            end

            @info "Repository updated successfully"
            return true
        end

    catch e
        @error "Error updating repository: $e"
        return false
    end
end

"""
    try_fallback_fetch(manager::RepositoryManager)

Try fetching from fallback remotes.
"""
function try_fallback_fetch(manager::RepositoryManager)
    for (i, fallback_url) in enumerate(manager.config.fallback_urls)
        remote_name = "fallback_$i"

        try
            # Add fallback remote if not exists
            remotes = read(`git remote`, String)
            if !occursin(remote_name, remotes)
                run(`git remote add $remote_name $fallback_url`)
            end

            # Try to fetch from fallback
            run(`git fetch $remote_name`)
            @info "Successfully fetched from fallback: $fallback_url"
            return true

        catch e
            @warn "Failed to fetch from fallback $fallback_url: $e"
            continue
        end
    end

    return false
end

"""
    verify_repository_integrity(manager::RepositoryManager)

Verify repository integrity and structure.
"""
function verify_repository_integrity(manager::RepositoryManager)
    @info "Verifying repository integrity..."

    integrity_checks = [
        check_git_repository,
        check_required_files,
        check_matlab_files,
        check_directory_structure
    ]

    all_passed = true

    for check in integrity_checks
        try
            result = check(manager.config.local_path)
            if !result
                all_passed = false
                @warn "Integrity check failed: $(nameof(check))"
            end
        catch e
            all_passed = false
            @warn "Integrity check error in $(nameof(check)): $e"
        end
    end

    manager.integrity_verified = all_passed
    manager.status["integrity_ok"] = all_passed

    if all_passed
        @info "Repository integrity verification passed"
    else
        @warn "Repository integrity verification failed"
    end

    return all_passed
end

"""
    check_git_repository(repo_path::String)

Check if directory is a valid git repository.
"""
function check_git_repository(repo_path::String)
    git_dir = joinpath(repo_path, ".git")

    if !isdir(git_dir)
        return false
    end

    try
        cd(repo_path) do
            # Check if it's a valid git repo
            run(`git status`)

            # Check if it has the correct remote
            remotes = read(`git remote -v`, String)
            return occursin("Helmholtz-adjoint-solver", remotes)
        end
    catch
        return false
    end
end

"""
    check_required_files(repo_path::String)

Check for required files in repository.
"""
function check_required_files(repo_path::String)
    required_files = ["README.md", "LICENSE"]

    for file in required_files
        found = false
        for (root, dirs, files) in walkdir(repo_path)
            if file in files
                found = true
                break
            end
        end
        if !found
            @warn "Required file not found: $file"
            return false
        end
    end

    return true
end

"""
    check_matlab_files(repo_path::String)

Check for critical MATLAB files.
"""
function check_matlab_files(repo_path::String)
    critical_matlab_files = ["ConvergentBornSolver.m"]

    for file in critical_matlab_files
        found = false
        for (root, dirs, files) in walkdir(repo_path)
            if file in files
                found = true
                break
            end
        end
        if !found
            @warn "Critical MATLAB file not found: $file"
            return false
        end
    end

    return true
end

"""
    check_directory_structure(repo_path::String)

Check expected directory structure.
"""
function check_directory_structure(repo_path::String)
    expected_dirs = ["src", "example"]

    for dir in expected_dirs
        found = false
        for (root, dirs, files) in walkdir(repo_path)
            if dir in dirs
                found = true
                break
            end
        end
        if !found
            @info "Expected directory not found (may be ok): $dir"
        end
    end

    return true  # Non-critical check
end

"""
    pin_to_commit(manager::RepositoryManager, commit_hash::String)

Pin repository to a specific commit.
"""
function pin_to_commit(manager::RepositoryManager, commit_hash::String)
    try
        cd(manager.config.local_path) do
            @info "Pinning repository to commit: $commit_hash"
            run(`git checkout $commit_hash`)
            @info "Repository pinned successfully"
        end
    catch e
        @error "Failed to pin to commit $commit_hash: $e"
        throw(e)
    end
end

"""
    check_offline_repository(manager::RepositoryManager)

Check if offline repository cache is available.
"""
function check_offline_repository(manager::RepositoryManager)
    cache_path = joinpath(manager.config.cache_path, "repository.tar.gz")
    return isfile(cache_path)
end

"""
    setup_offline_repository(manager::RepositoryManager)

Set up repository from offline cache.
"""
function setup_offline_repository(manager::RepositoryManager)
    cache_path = joinpath(manager.config.cache_path, "repository.tar.gz")

    if !isfile(cache_path)
        @error "Offline cache not found at $cache_path"
        return false
    end

    try
        @info "Extracting repository from offline cache..."

        # Create parent directory if needed
        mkpath(dirname(manager.config.local_path))

        # Extract cache
        run(`tar -xzf $cache_path -C $(dirname(manager.config.local_path))`)

        @info "Repository extracted from offline cache"
        return true

    catch e
        @error "Failed to extract offline repository: $e"
        return false
    end
end

"""
    create_offline_cache(manager::RepositoryManager)

Create offline cache of current repository.
"""
function create_offline_cache(manager::RepositoryManager)
    if !isdir(manager.config.local_path)
        @error "Cannot create cache - repository not found"
        return false
    end

    try
        @info "Creating offline repository cache..."

        # Create cache directory
        mkpath(manager.config.cache_path)

        cache_path = joinpath(manager.config.cache_path, "repository.tar.gz")

        # Create compressed archive
        run(`tar -czf $cache_path -C $(dirname(manager.config.local_path)) $(basename(manager.config.local_path))`)

        # Save metadata
        metadata = Dict(
            "created" => string(now()),
            "repository_url" => manager.config.primary_url,
            "branch" => manager.config.branch,
            "commit_hash" => get_current_commit(manager)
        )

        metadata_path = joinpath(manager.config.cache_path, "metadata.json")
        open(metadata_path, "w") do f
            JSON3.pretty(f, metadata)
        end

        @info "Offline cache created successfully"
        manager.offline_available = true
        return true

    catch e
        @error "Failed to create offline cache: $e"
        return false
    end
end

"""
    get_current_commit(manager::RepositoryManager)

Get current commit hash of repository.
"""
function get_current_commit(manager::RepositoryManager)
    if !isdir(manager.config.local_path)
        return "unknown"
    end

    try
        cd(manager.config.local_path) do
            return strip(read(`git rev-parse HEAD`, String))
        end
    catch
        return "unknown"
    end
end

"""
    get_repository_info(manager::RepositoryManager)

Get comprehensive repository information.
"""
function get_repository_info(manager::RepositoryManager)
    info = Dict{String, Any}()

    info["config"] = Dict(
        "primary_url" => manager.config.primary_url,
        "local_path" => manager.config.local_path,
        "branch" => manager.config.branch,
        "auto_update" => manager.config.auto_update,
        "offline_mode" => manager.config.offline_mode
    )

    info["status"] = copy(manager.status)
    info["last_update"] = manager.last_update
    info["integrity_verified"] = manager.integrity_verified
    info["fallback_active"] = manager.fallback_active
    info["offline_available"] = manager.offline_available

    if isdir(manager.config.local_path)
        try
            cd(manager.config.local_path) do
                info["current_commit"] = strip(read(`git rev-parse HEAD`, String))
                info["current_branch"] = strip(read(`git branch --show-current`, String))
                info["remote_url"] = strip(read(`git remote get-url origin`, String))

                # Check for pending updates
                try
                    run(`git fetch origin --dry-run`)
                    local_commit = strip(read(`git rev-parse HEAD`, String))
                    remote_commit = strip(read(`git rev-parse origin/$(manager.config.branch)`, String))
                    info["updates_available"] = local_commit != remote_commit
                catch
                    info["updates_available"] = "unknown"
                end
            end
        catch e
            info["git_error"] = string(e)
        end
    end

    return info
end

"""
    setup_fallback_mirrors(manager::RepositoryManager, mirror_urls::Vector{String})

Set up fallback mirror repositories.
"""
function setup_fallback_mirrors(manager::RepositoryManager, mirror_urls::Vector{String})
    @info "Setting up $(length(mirror_urls)) fallback mirrors"

    manager.config = RepositoryConfig(
        primary_url = manager.config.primary_url,
        fallback_urls = mirror_urls,
        local_path = manager.config.local_path,
        branch = manager.config.branch,
        commit_hash = manager.config.commit_hash,
        auto_update = manager.config.auto_update,
        verify_integrity = manager.config.verify_integrity,
        offline_mode = manager.config.offline_mode,
        cache_path = manager.config.cache_path
    )

    @info "Fallback mirrors configured successfully"
end

"""
    enable_offline_mode(manager::RepositoryManager)

Enable offline mode with cache creation.
"""
function enable_offline_mode(manager::RepositoryManager)
    @info "Enabling offline mode..."

    # Create cache if repository exists
    if isdir(manager.config.local_path)
        success = create_offline_cache(manager)
        if !success
            @warn "Failed to create offline cache"
            return false
        end
    end

    # Update configuration
    manager.config = RepositoryConfig(
        primary_url = manager.config.primary_url,
        fallback_urls = manager.config.fallback_urls,
        local_path = manager.config.local_path,
        branch = manager.config.branch,
        commit_hash = manager.config.commit_hash,
        auto_update = manager.config.auto_update,
        verify_integrity = manager.config.verify_integrity,
        offline_mode = true,
        cache_path = manager.config.cache_path
    )

    @info "Offline mode enabled"
    return true
end

"""
    sync_repository(manager::RepositoryManager; force_update=false)

Synchronize repository with comprehensive error handling.
"""
function sync_repository(manager::RepositoryManager; force_update = false)
    @info "Synchronizing repository..."

    # Check if repository exists
    if !isdir(manager.config.local_path)
        @info "Repository not found, cloning..."
        return clone_repository(manager)
    end

    # Update if auto-update is enabled or forced
    if manager.config.auto_update || force_update
        @info "Updating repository..."
        success = update_repository(manager)

        if !success
            @warn "Update failed, repository may be out of date"
            return false
        end
    end

    # Verify integrity
    if manager.config.verify_integrity
        integrity_ok = verify_repository_integrity(manager)
        if !integrity_ok
            @warn "Repository integrity check failed"
            return false
        end
    end

    @info "Repository synchronization completed successfully"
    return true
end

"""
    validate_repository_structure(repo_path::String)

Validate the complete repository structure for compatibility.
"""
function validate_repository_structure(repo_path::String)
    validation_results = Dict{String, Any}()

    # Check basic structure
    validation_results["is_directory"] = isdir(repo_path)
    validation_results["is_git_repo"] = check_git_repository(repo_path)
    validation_results["has_required_files"] = check_required_files(repo_path)
    validation_results["has_matlab_files"] = check_matlab_files(repo_path)
    validation_results["has_expected_structure"] = check_directory_structure(repo_path)

    # Check file counts
    if isdir(repo_path)
        matlab_files = []
        example_files = []

        for (root, dirs, files) in walkdir(repo_path)
            for file in files
                if endswith(file, ".m")
                    push!(matlab_files, relpath(joinpath(root, file), repo_path))
                elseif occursin("example", lowercase(root))
                    push!(example_files, relpath(joinpath(root, file), repo_path))
                end
            end
        end

        validation_results["matlab_file_count"] = length(matlab_files)
        validation_results["example_file_count"] = length(example_files)
        validation_results["matlab_files"] = matlab_files
        validation_results["example_files"] = example_files
    end

    # Overall validation
    critical_checks = [
        "is_directory", "is_git_repo", "has_required_files", "has_matlab_files"]
    validation_results["validation_passed"] = all(get(validation_results, check, false)
    for check in critical_checks)

    return validation_results
end

end # module RepositoryManager
