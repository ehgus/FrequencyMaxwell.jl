"""
MATLAB Integration Module for Cross-Validation Framework

Simplified MATLAB integration that just uses MATLAB.jl directly 
without environment variable modifications or complex detection logic.

Key features:
- Direct MATLAB.jl import without modifications
- Simple session management
- Basic error handling
"""

module MatlabIntegration

using Dates
using Logging
using UUIDs
using MATLAB

export initialize_matlab_session, cleanup_matlab_session
export MatlabSessionManager, is_matlab_available

"""
MATLAB.jl availability - true since this module requires MATLAB.jl
"""
const MATLAB_AVAILABLE = true

"""
Simple MATLAB session manager
"""
mutable struct MatlabSessionManager
    session::Union{Nothing, Any}  # MSession from MATLAB.jl
    session_id::String
    created_at::DateTime
    last_used::DateTime
    is_active::Bool
    
    function MatlabSessionManager()
        new(nothing, string(uuid4()), now(), now(), false)
    end
end

"""
    is_matlab_available()

Check if MATLAB.jl is available and working.
"""
function is_matlab_available()
    return MATLAB_AVAILABLE
end

"""
    initialize_matlab_session()

Initialize a simple MATLAB session using MATLAB.jl.
Returns a MatlabSessionManager with an active session.
"""
function initialize_matlab_session()
    if !MATLAB_AVAILABLE
        error("MATLAB.jl is not available. Please install it with: using Pkg; Pkg.add(\"MATLAB\")")
    end
    
    @info "Initializing MATLAB session..."
    
    try
        # Import MATLAB functions (already imported globally)
        session = MATLAB.MSession()
        
        # Test basic functionality
        MATLAB.eval_string(session, "fprintf('MATLAB session initialized successfully\\n')")
        
        # Create session manager
        manager = MatlabSessionManager()
        manager.session = session
        manager.is_active = true
        manager.last_used = now()
        
        @info "MATLAB session initialized successfully"
        return manager
        
    catch e
        @error "Failed to initialize MATLAB session" exception=e
        rethrow(e)
    end
end

"""
    cleanup_matlab_session(manager::MatlabSessionManager)

Clean up MATLAB session and resources.
"""
function cleanup_matlab_session(manager::MatlabSessionManager)
    if manager.is_active && manager.session !== nothing
        try
            # Close MATLAB session
            MATLAB.close(manager.session)
            manager.is_active = false
            manager.session = nothing
            @info "MATLAB session cleaned up successfully"
        catch e
            @warn "Error during MATLAB session cleanup" exception=e
        end
    end
end

end # module MatlabIntegration