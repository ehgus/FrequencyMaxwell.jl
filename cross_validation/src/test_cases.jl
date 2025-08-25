"""
Predefined Test Cases for FrequencyMaxwell Cross-Validation Framework

This file contains predefined test case functions that can be used
with the cross-validation framework. Each function returns a TestCase
object with appropriate parameters and tolerances.

Usage:
    using CrossValidation
    cv = CrossValidationFramework()
    add_test_case!(cv, metalens_forward_test())
    add_test_case!(cv, grating_adjoint_test())
"""

"""
    metalens_forward_test(; grid_size=[64, 64], wavelength=532e-9, focal_length=50e-6)

Create a forward solver test case for metalens validation.
"""
function metalens_forward_test(; 
    grid_size=[64, 64], 
    wavelength=532e-9, 
    focal_length=50e-6,
    tolerance_level="standard"
)
    
    # Set tolerances based on level
    if tolerance_level == "strict"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-8,
            "energy_relative_error" => 1e-10,
            "convergence_threshold" => 1e-12,
            "focusing_efficiency_relative_error" => 1e-6
        )
    elseif tolerance_level == "relaxed"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-3,
            "energy_relative_error" => 1e-5,
            "convergence_threshold" => 1e-6,
            "focusing_efficiency_relative_error" => 1e-2
        )
    else  # standard
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-6,
            "energy_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "focusing_efficiency_relative_error" => 1e-4
        )
    end
    
    parameters = Dict{String, Any}(
        "wavelength" => wavelength,
        "focal_length" => focal_length,
        "grid_size" => grid_size,
        "test_name" => "metalens_forward"
    )
    
    description = "Forward solver validation for metalens focusing with " *
                 "$(grid_size[1])×$(grid_size[2]) grid at $(wavelength*1e9) nm"
    
    return TestCase(
        "metalens_forward",
        "forward",
        parameters,
        tolerances,
        description
    )
end

"""
    grating_forward_test(; grid_size=[128, 32], wavelength=532e-9, grating_period=2e-6)

Create a forward solver test case for diffraction grating validation.
"""
function grating_forward_test(; 
    grid_size=[128, 32], 
    wavelength=532e-9, 
    grating_period=2e-6,
    tolerance_level="standard"
)
    
    if tolerance_level == "strict"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-8,
            "energy_relative_error" => 1e-10,
            "convergence_threshold" => 1e-12,
            "efficiency_relative_error" => 1e-6
        )
    elseif tolerance_level == "relaxed"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-3,
            "energy_relative_error" => 1e-5,
            "convergence_threshold" => 1e-6,
            "efficiency_relative_error" => 1e-2
        )
    else  # standard
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-6,
            "energy_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "efficiency_relative_error" => 1e-4
        )
    end
    
    parameters = Dict{String, Any}(
        "wavelength" => wavelength,
        "grating_period" => grating_period,
        "grid_size" => grid_size,
        "test_name" => "grating_forward"
    )
    
    description = "Forward solver validation for diffraction grating with " *
                 "period $(grating_period*1e6) μm at $(wavelength*1e9) nm"
    
    return TestCase(
        "grating_forward",
        "forward",
        parameters,
        tolerances,
        description
    )
end

"""
    two_beam_forward_test(; grid_size=[64, 64], wavelength=532e-9, beam_angle=15.0)

Create a forward solver test case for two-beam interference validation.
"""
function two_beam_forward_test(; 
    grid_size=[64, 64], 
    wavelength=532e-9, 
    beam_angle=15.0,
    tolerance_level="standard"
)
    
    if tolerance_level == "strict"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-8,
            "energy_relative_error" => 1e-10,
            "convergence_threshold" => 1e-12,
            "contrast_relative_error" => 1e-6
        )
    elseif tolerance_level == "relaxed"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-3,
            "energy_relative_error" => 1e-5,
            "convergence_threshold" => 1e-6,
            "contrast_relative_error" => 1e-2
        )
    else  # standard
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-6,
            "energy_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "contrast_relative_error" => 1e-4
        )
    end
    
    parameters = Dict{String, Any}(
        "wavelength" => wavelength,
        "beam_angle" => beam_angle,
        "grid_size" => grid_size,
        "test_name" => "two_beam_forward"
    )
    
    description = "Two-beam interference forward validation with " *
                 "$(beam_angle)° beam angle at $(wavelength*1e9) nm"
    
    return TestCase(
        "two_beam_forward",
        "forward",
        parameters,
        tolerances,
        description
    )
end

"""
    sio2_sphere_forward_test(; grid_size=[64, 64], wavelength=532e-9, sphere_radius=2.5e-6)

Create a forward solver test case for SiO2 sphere scattering validation.
"""
function sio2_sphere_forward_test(; 
    grid_size=[64, 64], 
    wavelength=532e-9, 
    sphere_radius=2.5e-6,
    tolerance_level="standard"
)
    
    if tolerance_level == "strict"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-8,
            "energy_relative_error" => 1e-10,
            "convergence_threshold" => 1e-12,
            "scattering_cross_section_relative_error" => 1e-6
        )
    elseif tolerance_level == "relaxed"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-3,
            "energy_relative_error" => 1e-5,
            "convergence_threshold" => 1e-6,
            "scattering_cross_section_relative_error" => 1e-2
        )
    else  # standard
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-6,
            "energy_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "scattering_cross_section_relative_error" => 1e-4
        )
    end
    
    parameters = Dict{String, Any}(
        "wavelength" => wavelength,
        "sphere_radius" => sphere_radius,
        "grid_size" => grid_size,
        "test_name" => "sio2_sphere_forward"
    )
    
    description = "SiO2 sphere scattering validation with " *
                 "$(sphere_radius*1e6) μm radius at $(wavelength*1e9) nm"
    
    return TestCase(
        "sio2_sphere_forward",
        "forward",
        parameters,
        tolerances,
        description
    )
end

"""
    helical_metalens_forward_test(; grid_size=[64, 64], wavelength=532e-9, topological_charge=1, focal_length=50e-6)

Create a forward solver test case for helical metalens validation.
"""
function helical_metalens_forward_test(; 
    grid_size=[64, 64], 
    wavelength=532e-9, 
    topological_charge=1, 
    focal_length=50e-6,
    tolerance_level="standard"
)
    
    if tolerance_level == "strict"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-8,
            "energy_relative_error" => 1e-10,
            "convergence_threshold" => 1e-12,
            "measured_topological_charge_relative_error" => 1e-6
        )
    elseif tolerance_level == "relaxed"
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-3,
            "energy_relative_error" => 1e-5,
            "convergence_threshold" => 1e-6,
            "measured_topological_charge_relative_error" => 1e-1
        )
    else  # standard
        tolerances = Dict{String, Float64}(
            "field_relative_error" => 1e-6,
            "energy_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "measured_topological_charge_relative_error" => 1e-2
        )
    end
    
    parameters = Dict{String, Any}(
        "wavelength" => wavelength,
        "topological_charge" => topological_charge,
        "focal_length" => focal_length,
        "grid_size" => grid_size,
        "test_name" => "helical_metalens_forward"
    )
    
    description = "Helical metalens forward validation with charge " *
                 "$topological_charge and $(focal_length*1e6) μm focal length"
    
    return TestCase(
        "helical_metalens_forward",
        "forward",
        parameters,
        tolerances,
        description
    )
end

# Adjoint Test Cases

"""
    metalens_adjoint_test(; grid_size=[32, 32], wavelength=532e-9, target_focal_length=50e-6, max_iterations=10)

Create an adjoint solver test case for metalens optimization.
"""
function metalens_adjoint_test(; 
    grid_size=[32, 32], 
    wavelength=532e-9, 
    target_focal_length=50e-6,
    max_iterations=10,
    tolerance_level="standard"
)
    
    if tolerance_level == "strict"
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "focusing_efficiency_relative_error" => 1e-6,
            "final_residual" => 1e-10
        )
    elseif tolerance_level == "relaxed"
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-2,
            "convergence_threshold" => 1e-4,
            "focusing_efficiency_relative_error" => 1e-1,
            "final_residual" => 1e-4
        )
    else  # standard
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-5,
            "convergence_threshold" => 1e-8,
            "focusing_efficiency_relative_error" => 1e-3,
            "final_residual" => 1e-6
        )
    end
    
    parameters = Dict{String, Any}(
        "wavelength" => wavelength,
        "target_focal_length" => target_focal_length,
        "grid_size" => grid_size,
        "max_iterations" => max_iterations,
        "test_name" => "metalens_adjoint"
    )
    
    description = "Metalens adjoint optimization with $(grid_size[1])×$(grid_size[2]) grid, " *
                 "$max_iterations iterations"
    
    return TestCase(
        "metalens_adjoint",
        "adjoint",
        parameters,
        tolerances,
        description
    )
end

"""
    grating_adjoint_test(; grid_size=[64, 16], wavelength=532e-9, target_efficiency=0.8, max_iterations=15)

Create an adjoint solver test case for grating efficiency optimization.
"""
function grating_adjoint_test(; 
    grid_size=[64, 16], 
    wavelength=532e-9, 
    target_efficiency=0.8,
    max_iterations=15,
    tolerance_level="standard"
)
    
    if tolerance_level == "strict"
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "final_efficiency_relative_error" => 1e-6,
            "final_residual" => 1e-10
        )
    elseif tolerance_level == "relaxed"
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-2,
            "convergence_threshold" => 1e-4,
            "final_efficiency_relative_error" => 1e-1,
            "final_residual" => 1e-4
        )
    else  # standard
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-5,
            "convergence_threshold" => 1e-8,
            "final_efficiency_relative_error" => 1e-3,
            "final_residual" => 1e-6
        )
    end
    
    parameters = Dict{String, Any}(
        "wavelength" => wavelength,
        "target_efficiency" => target_efficiency,
        "grid_size" => grid_size,
        "max_iterations" => max_iterations,
        "test_name" => "grating_adjoint"
    )
    
    description = "Grating adjoint optimization targeting " *
                 "$(target_efficiency*100)% efficiency with $max_iterations iterations"
    
    return TestCase(
        "grating_adjoint",
        "adjoint",
        parameters,
        tolerances,
        description
    )
end

"""
    double_helix_adjoint_test(; grid_size=[48, 48], wavelength=532e-9, helix_separation=5e-6, max_iterations=8)

Create an adjoint solver test case for double helix PSF optimization.
"""
function double_helix_adjoint_test(; 
    grid_size=[48, 48], 
    wavelength=532e-9, 
    helix_separation=5e-6,
    max_iterations=8,
    tolerance_level="standard"
)
    
    if tolerance_level == "strict"
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "contrast_ratio_relative_error" => 1e-6,
            "final_residual" => 1e-10
        )
    elseif tolerance_level == "relaxed"
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-2,
            "convergence_threshold" => 1e-4,
            "contrast_ratio_relative_error" => 1e-1,
            "final_residual" => 1e-4
        )
    else  # standard
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-5,
            "convergence_threshold" => 1e-8,
            "contrast_ratio_relative_error" => 1e-3,
            "final_residual" => 1e-6
        )
    end
    
    parameters = Dict{String, Any}(
        "wavelength" => wavelength,
        "helix_separation" => helix_separation,
        "grid_size" => grid_size,
        "max_iterations" => max_iterations,
        "test_name" => "double_helix_adjoint"
    )
    
    description = "Double helix PSF optimization with " *
                 "$(helix_separation*1e6) μm separation and $max_iterations iterations"
    
    return TestCase(
        "double_helix_adjoint",
        "adjoint",
        parameters,
        tolerances,
        description
    )
end

"""
    single_helix_adjoint_test(; grid_size=[40, 40], wavelength=532e-9, topological_charge=2, max_iterations=12)

Create an adjoint solver test case for single helix PSF optimization.
"""
function single_helix_adjoint_test(; 
    grid_size=[40, 40], 
    wavelength=532e-9, 
    topological_charge=2,
    max_iterations=12,
    tolerance_level="standard"
)
    
    if tolerance_level == "strict"
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-8,
            "convergence_threshold" => 1e-10,
            "rotation_response_relative_error" => 1e-6,
            "final_residual" => 1e-10
        )
    elseif tolerance_level == "relaxed"
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-2,
            "convergence_threshold" => 1e-4,
            "rotation_response_relative_error" => 1e-1,
            "final_residual" => 1e-4
        )
    else  # standard
        tolerances = Dict{String, Float64}(
            "gradient_relative_error" => 1e-5,
            "convergence_threshold" => 1e-8,
            "rotation_response_relative_error" => 1e-3,
            "final_residual" => 1e-6
        )
    end
    
    parameters = Dict{String, Any}(
        "wavelength" => wavelength,
        "topological_charge" => topological_charge,
        "grid_size" => grid_size,
        "max_iterations" => max_iterations,
        "test_name" => "single_helix_adjoint"
    )
    
    description = "Single helix PSF optimization with charge " *
                 "$topological_charge and $max_iterations iterations"
    
    return TestCase(
        "single_helix_adjoint",
        "adjoint",
        parameters,
        tolerances,
        description
    )
end

# Utility functions for creating test suites

"""
    create_forward_test_suite(; tolerance_level="standard")

Create a comprehensive suite of forward solver test cases.
"""
function create_forward_test_suite(; tolerance_level="standard")
    return [
        metalens_forward_test(tolerance_level=tolerance_level),
        grating_forward_test(tolerance_level=tolerance_level),
        two_beam_forward_test(tolerance_level=tolerance_level),
        sio2_sphere_forward_test(tolerance_level=tolerance_level),
        helical_metalens_forward_test(tolerance_level=tolerance_level)
    ]
end

"""
    create_adjoint_test_suite(; tolerance_level="standard")

Create a comprehensive suite of adjoint solver test cases.
"""
function create_adjoint_test_suite(; tolerance_level="standard")
    return [
        metalens_adjoint_test(tolerance_level=tolerance_level),
        grating_adjoint_test(tolerance_level=tolerance_level),
        double_helix_adjoint_test(tolerance_level=tolerance_level),
        single_helix_adjoint_test(tolerance_level=tolerance_level)
    ]
end

"""
    create_full_test_suite(; tolerance_level="standard")

Create a comprehensive suite including both forward and adjoint test cases.
"""
function create_full_test_suite(; tolerance_level="standard")
    forward_tests = create_forward_test_suite(tolerance_level=tolerance_level)
    adjoint_tests = create_adjoint_test_suite(tolerance_level=tolerance_level)
    return [forward_tests; adjoint_tests]
end

"""
    create_quick_test_suite(; tolerance_level="relaxed")

Create a quick test suite with smaller problem sizes for development and debugging.
"""
function create_quick_test_suite(; tolerance_level="relaxed")
    return [
        metalens_forward_test(grid_size=[16, 16], tolerance_level=tolerance_level),
        grating_forward_test(grid_size=[32, 8], tolerance_level=tolerance_level),
        metalens_adjoint_test(grid_size=[16, 16], max_iterations=3, tolerance_level=tolerance_level)
    ]
end

"""
    create_performance_test_suite(; tolerance_level="standard")

Create a performance test suite with larger problem sizes.
"""
function create_performance_test_suite(; tolerance_level="standard")
    return [
        metalens_forward_test(grid_size=[128, 128], tolerance_level=tolerance_level),
        grating_forward_test(grid_size=[256, 64], tolerance_level=tolerance_level),
        helical_metalens_forward_test(grid_size=[128, 128], tolerance_level=tolerance_level),
        metalens_adjoint_test(grid_size=[64, 64], max_iterations=20, tolerance_level=tolerance_level),
        grating_adjoint_test(grid_size=[128, 32], max_iterations=25, tolerance_level=tolerance_level)
    ]
end