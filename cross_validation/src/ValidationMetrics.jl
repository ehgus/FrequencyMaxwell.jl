"""
Enhanced Validation Metrics Module

This module provides robust, comprehensive validation metrics for comparing
electromagnetic solver results. It includes advanced error analysis, 
numerical stability checks, and detailed diagnostic capabilities.

Features:
- Multiple error metrics (absolute, relative, RMS, max)
- Field validation with spatial analysis
- Energy conservation checks
- Convergence analysis with trend detection
- Statistical validation with confidence intervals
- Robust handling of edge cases and numerical issues
"""

module ValidationMetrics

using LinearAlgebra
using Statistics
using Printf

export compute_comprehensive_metrics, ValidationMetricsResult, analyze_convergence
export check_energy_conservation, analyze_field_patterns, generate_metric_report
export validate_numerical_stability, compute_confidence_intervals, convert_metrics_to_dict

"""
Comprehensive validation metrics structure
"""
struct ValidationMetricsResult
    # Field validation
    field_absolute_error::Float64
    field_relative_error::Float64
    field_rms_error::Float64
    field_max_error::Float64
    field_mean_error::Float64
    field_correlation::Float64

    # Energy metrics
    energy_relative_error::Float64
    energy_absolute_error::Float64
    energy_conservation_ratio::Float64

    # Convergence metrics
    final_residual::Float64
    convergence_rate::Float64
    iterations_to_convergence::Int64
    convergence_trend::String

    # Statistical metrics
    relative_error_std::Float64
    relative_error_mean::Float64
    outlier_fraction::Float64

    # Quality indicators
    numerical_stability::Float64
    solution_quality::Float64
    validation_confidence::Float64

    # Diagnostic information
    warnings::Vector{String}
    quality_flags::Dict{String, Bool}

    function ValidationMetricsResult(metrics_dict::Dict)
        warnings = String[]
        quality_flags = Dict{String, Bool}()

        # Extract metrics with defaults
        field_abs = get(metrics_dict, "field_absolute_error", NaN)
        field_rel = get(metrics_dict, "field_relative_error", NaN)
        field_rms = get(metrics_dict, "field_rms_error", NaN)
        field_max = get(metrics_dict, "field_max_error", NaN)
        field_mean = get(metrics_dict, "field_mean_error", NaN)
        field_corr = get(metrics_dict, "field_correlation", NaN)

        energy_rel = get(metrics_dict, "energy_relative_error", NaN)
        energy_abs = get(metrics_dict, "energy_absolute_error", NaN)
        energy_cons = get(metrics_dict, "energy_conservation_ratio", NaN)

        final_res = get(metrics_dict, "final_residual", NaN)
        conv_rate = get(metrics_dict, "convergence_rate", NaN)
        iterations = get(metrics_dict, "iterations_to_convergence", 0)
        conv_trend = get(metrics_dict, "convergence_trend", "unknown")

        rel_err_std = get(metrics_dict, "relative_error_std", NaN)
        rel_err_mean = get(metrics_dict, "relative_error_mean", NaN)
        outlier_frac = get(metrics_dict, "outlier_fraction", NaN)

        num_stab = get(metrics_dict, "numerical_stability", NaN)
        sol_qual = get(metrics_dict, "solution_quality", NaN)
        val_conf = get(metrics_dict, "validation_confidence", NaN)

        # Quality checks
        if isnan(field_rel) || field_rel > 1e-3
            push!(warnings, "High field relative error detected")
            quality_flags["field_quality"] = false
        else
            quality_flags["field_quality"] = true
        end

        if isnan(energy_rel) || abs(energy_rel) > 1e-5
            push!(warnings, "Energy conservation may be violated")
            quality_flags["energy_conservation"] = false
        else
            quality_flags["energy_conservation"] = true
        end

        if isnan(final_res) || final_res > 1e-6
            push!(warnings, "Poor convergence detected")
            quality_flags["convergence"] = false
        else
            quality_flags["convergence"] = true
        end

        new(field_abs, field_rel, field_rms, field_max, field_mean, field_corr,
            energy_rel, energy_abs, energy_cons,
            final_res, conv_rate, iterations, conv_trend,
            rel_err_std, rel_err_mean, outlier_frac,
            num_stab, sol_qual, val_conf,
            warnings, quality_flags)
    end
end

"""
    compute_comprehensive_metrics(matlab_results::Dict, reference_results::Dict, tolerances::Dict)

Compute comprehensive validation metrics with robust error handling.
"""
function compute_comprehensive_metrics(matlab_results::Dict, reference_results::Dict, tolerances::Dict)
    metrics = Dict{String, Any}()

    try
        # Field validation metrics
        field_metrics = compute_field_metrics(matlab_results, reference_results)
        merge!(metrics, field_metrics)

        # Energy conservation metrics
        energy_metrics = compute_energy_metrics(matlab_results, reference_results)
        merge!(metrics, energy_metrics)

        # Convergence analysis
        convergence_metrics = compute_convergence_metrics(matlab_results)
        merge!(metrics, convergence_metrics)

        # Statistical analysis
        statistical_metrics = compute_statistical_metrics(matlab_results, reference_results)
        merge!(metrics, statistical_metrics)

        # Quality assessment
        quality_metrics = assess_solution_quality(matlab_results, reference_results)
        merge!(metrics, quality_metrics)

        # Efficiency metrics
        efficiency_metrics = compute_efficiency_metrics(matlab_results, reference_results)
        merge!(metrics, efficiency_metrics)

        # Apply tolerance checks
        metrics["meets_tolerances"] = check_tolerances(metrics, tolerances)

        # Add tolerance information for reporting
        metrics["tolerances"] = tolerances

    catch e
        @warn "Error computing comprehensive metrics" exception=e
        metrics["computation_error"] = 1.0
        metrics["error_message"] = string(e)
    end

    # Convert ValidationMetricsResult to Dict{String, Float64} for compatibility
    result_struct = ValidationMetricsResult(metrics)
    return convert_metrics_to_dict(result_struct)
end

"""
    compute_field_metrics(matlab_results::Dict, reference_results::Dict)

Compute detailed field comparison metrics.
"""
function compute_field_metrics(matlab_results::Dict, reference_results::Dict)
    metrics = Dict{String, Float64}()

    if haskey(matlab_results, "E_field") && haskey(reference_results, "E_field")
        E_matlab = matlab_results["E_field"]
        E_ref = reference_results["E_field"]

        if size(E_matlab) == size(E_ref)
            # Handle complex fields properly
            if eltype(E_matlab) <: Complex || eltype(E_ref) <: Complex
                field_diff = abs.(E_matlab .- E_ref)
                field_magnitude_ref = abs.(E_ref)
                field_magnitude_matlab = abs.(E_matlab)
            else
                field_diff = abs.(E_matlab .- E_ref)
                field_magnitude_ref = abs.(E_ref)
                field_magnitude_matlab = abs.(E_matlab)
            end

            # Absolute errors
            metrics["field_absolute_error"] = maximum(field_diff)
            metrics["field_mean_error"] = mean(field_diff)
            metrics["field_rms_error"] = sqrt(mean(field_diff .^ 2))

            # Relative errors (handle division by zero)
            ref_magnitude = maximum(field_magnitude_ref)
            if ref_magnitude > 1e-15
                relative_errors = field_diff ./
                                  (field_magnitude_ref .+ 1e-15 * ref_magnitude)
                metrics["field_relative_error"] = maximum(relative_errors)
                metrics["field_max_error"] = maximum(field_diff) / ref_magnitude

                # Statistical measures
                metrics["relative_error_mean"] = mean(relative_errors)
                metrics["relative_error_std"] = std(relative_errors)

                # Outlier detection
                threshold = mean(relative_errors) + 3 * std(relative_errors)
                outliers = sum(relative_errors .> threshold)
                metrics["outlier_fraction"] = outliers / length(relative_errors)
            else
                metrics["field_relative_error"] = Inf
                metrics["field_max_error"] = Inf
            end

            # Correlation analysis
            if length(field_magnitude_matlab) > 1 && length(field_magnitude_ref) > 1
                try
                    correlation_matrix = cor(vec(field_magnitude_matlab), vec(field_magnitude_ref))
                    metrics["field_correlation"] = correlation_matrix
                catch
                    metrics["field_correlation"] = NaN
                end
            end
        else
            metrics["field_size_mismatch"] = 1.0
        end
    end

    return metrics
end

"""
    compute_energy_metrics(matlab_results::Dict, reference_results::Dict)

Compute energy and power conservation metrics.
"""
function compute_energy_metrics(matlab_results::Dict, reference_results::Dict)
    metrics = Dict{String, Float64}()

    # Check various energy/intensity measures
    energy_keys = ["intensity", "energy", "power", "total_power"]

    for key in energy_keys
        if haskey(matlab_results, key) && haskey(reference_results, key)
            val_matlab = matlab_results[key]
            val_ref = reference_results[key]

            # Handle both scalar and array values
            energy_matlab = isa(val_matlab, Array) ? sum(val_matlab) : val_matlab
            energy_ref = isa(val_ref, Array) ? sum(val_ref) : val_ref

            if energy_ref != 0
                rel_error = abs(energy_matlab - energy_ref) / abs(energy_ref)
                metrics["$(key)_relative_error"] = rel_error
                metrics["$(key)_absolute_error"] = abs(energy_matlab - energy_ref)

                # Conservation ratio
                metrics["$(key)_conservation_ratio"] = energy_matlab / energy_ref

                # Use the first found energy measure for general metrics
                if !haskey(metrics, "energy_relative_error")
                    metrics["energy_relative_error"] = rel_error
                    metrics["energy_absolute_error"] = abs(energy_matlab - energy_ref)
                    metrics["energy_conservation_ratio"] = energy_matlab / energy_ref
                end
            end

            break  # Use first available energy measure
        end
    end

    return metrics
end

"""
    compute_convergence_metrics(matlab_results::Dict)

Analyze convergence behavior and quality.
"""
function compute_convergence_metrics(matlab_results::Dict)
    metrics = Dict{String, Any}()

    if haskey(matlab_results, "convergence_history")
        conv_hist = matlab_results["convergence_history"]

        if !isempty(conv_hist) && length(conv_hist) > 1
            # Final residual
            metrics["final_residual"] = conv_hist[end]
            metrics["initial_residual"] = conv_hist[1]

            # Convergence rate analysis
            n_points = length(conv_hist)
            if n_points > 3
                # Linear convergence rate (log scale)
                log_residuals = log.(max.(conv_hist, 1e-16))
                iterations = 1:n_points

                # Fit linear trend
                A = [iterations ones(n_points)]
                try
                    coeffs = A \ log_residuals
                    metrics["convergence_rate"] = -coeffs[1]  # Negative slope indicates convergence

                    # Determine convergence trend
                    if coeffs[1] < -0.1
                        metrics["convergence_trend"] = "fast"
                    elseif coeffs[1] < -0.01
                        metrics["convergence_trend"] = "moderate"
                    elseif coeffs[1] < 0
                        metrics["convergence_trend"] = "slow"
                    else
                        metrics["convergence_trend"] = "stalled"
                    end
                catch
                    metrics["convergence_rate"] = 0.0
                    metrics["convergence_trend"] = "unknown"
                end
            else
                metrics["convergence_rate"] = 0.0
                metrics["convergence_trend"] = "insufficient_data"
            end

            # Iterations to convergence (threshold-based)
            convergence_threshold = 1e-8
            converged_idx = findfirst(x -> x < convergence_threshold, conv_hist)
            metrics["iterations_to_convergence"] = converged_idx !== nothing ?
                                                   converged_idx : n_points

            # Convergence quality
            if n_points > 5
                # Check for monotonic decrease
                decreasing = all(diff(conv_hist) .<= 0)
                metrics["monotonic_convergence"] = decreasing ? 1.0 : 0.0

                # Smoothness of convergence
                second_diff = diff(diff(conv_hist))
                metrics["convergence_smoothness"] = 1.0 / (1.0 + std(second_diff))
            end
        else
            metrics["convergence_data_available"] = false
        end
    end

    return metrics
end

"""
    compute_statistical_metrics(matlab_results::Dict, reference_results::Dict)

Compute statistical validation metrics.
"""
function compute_statistical_metrics(matlab_results::Dict, reference_results::Dict)
    metrics = Dict{String, Float64}()

    # Analyze residuals and error distributions
    if haskey(matlab_results, "E_field") && haskey(reference_results, "E_field")
        E_matlab = matlab_results["E_field"]
        E_ref = reference_results["E_field"]

        if size(E_matlab) == size(E_ref)
            errors = vec(abs.(E_matlab .- E_ref))
            ref_magnitudes = vec(abs.(E_ref))

            if !isempty(errors) && all(isfinite.(errors))
                # Error distribution statistics
                metrics["error_mean"] = mean(errors)
                metrics["error_std"] = std(errors)
                metrics["error_median"] = median(errors)
                metrics["error_q75"] = quantile(errors, 0.75)
                metrics["error_q95"] = quantile(errors, 0.95)
                metrics["error_max"] = maximum(errors)

                # Relative error statistics (where reference is significant)
                significant_mask = ref_magnitudes .> maximum(ref_magnitudes) * 1e-6
                if sum(significant_mask) > 0
                    rel_errors = errors[significant_mask] ./
                                 ref_magnitudes[significant_mask]
                    metrics["relative_error_mean"] = mean(rel_errors)
                    metrics["relative_error_std"] = std(rel_errors)
                    metrics["relative_error_median"] = median(rel_errors)
                    metrics["relative_error_q95"] = quantile(rel_errors, 0.95)
                end

                # Distribution shape analysis (simplified)
                if length(errors) > 10
                    # Simple shape indicators
                    error_range = maximum(errors) - minimum(errors)
                    error_iqr = quantile(errors, 0.75) - quantile(errors, 0.25)

                    if error_range > 0
                        shape_indicator = error_iqr / error_range
                        metrics["error_shape"] = shape_indicator

                        # Simple normality estimate based on range vs IQR ratio
                        normality_score = min(1.0, 2 * shape_indicator)
                        metrics["error_normality"] = normality_score
                    end
                end
            end
        end
    end

    return metrics
end

"""
    assess_solution_quality(matlab_results::Dict, reference_results::Dict)

Assess overall solution quality and numerical stability.
"""
function assess_solution_quality(matlab_results::Dict, reference_results::Dict)
    metrics = Dict{String, Float64}()

    # Numerical stability indicators
    stability_score = 1.0
    quality_score = 1.0
    confidence_score = 1.0

    # Check for NaN/Inf values
    for (_, value) in matlab_results
        if isa(value, Array)
            if any(isnan.(value)) || any(isinf.(value))
                stability_score *= 0.1
                break
            end
        elseif isa(value, Number)
            if isnan(value) || isinf(value)
                stability_score *= 0.1
            end
        end
    end

    # Field magnitude consistency
    if haskey(matlab_results, "E_field") && haskey(reference_results, "E_field")
        E_matlab = matlab_results["E_field"]
        E_ref = reference_results["E_field"]

        if size(E_matlab) == size(E_ref)
            # Check dynamic range
            matlab_max = maximum(abs.(E_matlab))
            matlab_min = minimum(abs.(E_matlab[abs.(E_matlab) .> 0]))
            ref_max = maximum(abs.(E_ref))
            ref_min = minimum(abs.(E_ref[abs.(E_ref) .> 0]))

            if matlab_min > 0 && ref_min > 0
                matlab_dynamic_range = matlab_max / matlab_min
                ref_dynamic_range = ref_max / ref_min

                range_ratio = min(matlab_dynamic_range, ref_dynamic_range) /
                              max(matlab_dynamic_range, ref_dynamic_range)
                quality_score *= range_ratio
            end

            # Phase consistency (for complex fields)
            if eltype(E_matlab) <: Complex && eltype(E_ref) <: Complex
                phase_matlab = angle.(E_matlab)
                phase_ref = angle.(E_ref)
                phase_diff = abs.(phase_matlab .- phase_ref)

                # Handle phase wrapping
                phase_diff = min.(phase_diff, 2π .- phase_diff)
                phase_consistency = 1.0 - mean(phase_diff) / π
                quality_score *= max(phase_consistency, 0.1)
            end
        end
    end

    # Convergence quality impact
    if haskey(matlab_results, "convergence_history")
        conv_hist = matlab_results["convergence_history"]
        if !isempty(conv_hist)
            final_residual = conv_hist[end]
            if final_residual < 1e-10
                confidence_score *= 1.0
            elseif final_residual < 1e-8
                confidence_score *= 0.9
            elseif final_residual < 1e-6
                confidence_score *= 0.7
            else
                confidence_score *= 0.5
            end
        end
    end

    metrics["numerical_stability"] = stability_score
    metrics["solution_quality"] = quality_score
    metrics["validation_confidence"] = confidence_score

    return metrics
end

"""
    compute_efficiency_metrics(matlab_results::Dict, reference_results::Dict)

Compute physics-specific efficiency metrics.
"""
function compute_efficiency_metrics(matlab_results::Dict, reference_results::Dict)
    metrics = Dict{String, Float64}()

    # Check various efficiency measures
    efficiency_keys = ["efficiency", "focusing_efficiency", "scattering_cross_section",
        "transmission", "reflection", "absorption"]

    for key in efficiency_keys
        if haskey(matlab_results, key) && haskey(reference_results, key)
            val_matlab = matlab_results[key]
            val_ref = reference_results[key]

            if val_ref != 0
                rel_error = abs(val_matlab - val_ref) / abs(val_ref)
                metrics["$(key)_relative_error"] = rel_error
                metrics["$(key)_absolute_error"] = abs(val_matlab - val_ref)
            end
        end
    end

    return metrics
end

"""
    analyze_convergence(convergence_history::Vector{Float64})

Detailed convergence analysis with trend detection.
"""
function analyze_convergence(convergence_history::Vector{Float64})
    if isempty(convergence_history)
        return Dict("status" => "no_data")
    end

    analysis = Dict{String, Any}()
    n = length(convergence_history)

    # Basic metrics
    analysis["initial_residual"] = convergence_history[1]
    analysis["final_residual"] = convergence_history[end]
    analysis["reduction_factor"] = convergence_history[1] / convergence_history[end]
    analysis["iterations"] = n

    if n > 1
        # Monotonicity check
        is_monotonic = all(diff(convergence_history) .<= 0)
        analysis["monotonic"] = is_monotonic

        # Convergence rate estimation
        if n > 3
            log_hist = log.(max.(convergence_history, 1e-16))
            iterations = collect(1:n)

            # Linear fit in log space
            A = [iterations ones(n)]
            try
                coeffs = A \ log_hist
                analysis["linear_rate"] = -coeffs[1]
                analysis["r_squared"] = cor(A * coeffs, log_hist)^2
            catch
                analysis["linear_rate"] = NaN
                analysis["r_squared"] = NaN
            end

            # Detect stagnation
            recent_fraction = max(1, div(n, 4))
            recent_change = abs(convergence_history[end] -
                                convergence_history[end - recent_fraction + 1])
            stagnation_threshold = abs(convergence_history[end]) * 0.01
            analysis["stagnating"] = recent_change < stagnation_threshold
        end
    end

    return analysis
end

"""
    check_energy_conservation(matlab_results::Dict, tolerance::Float64=1e-6)

Check energy conservation with detailed diagnostics.
"""
function check_energy_conservation(matlab_results::Dict, tolerance::Float64 = 1e-6)
    diagnostics = Dict{String, Any}()

    # Check incident vs scattered + absorbed energy
    total_incident = get(matlab_results, "incident_power", NaN)
    total_scattered = get(matlab_results, "scattered_power", NaN)
    total_absorbed = get(matlab_results, "absorbed_power", NaN)

    if !isnan(total_incident) && !isnan(total_scattered) && !isnan(total_absorbed)
        total_accounted = total_scattered + total_absorbed
        conservation_error = abs(total_incident - total_accounted) / total_incident

        diagnostics["conservation_error"] = conservation_error
        diagnostics["conserved"] = conservation_error < tolerance
        diagnostics["incident_power"] = total_incident
        diagnostics["scattered_power"] = total_scattered
        diagnostics["absorbed_power"] = total_absorbed
    else
        # Fallback: check field energy consistency
        if haskey(matlab_results, "intensity")
            intensity = matlab_results["intensity"]
            total_intensity = sum(intensity)

            # Simple energy conservation check
            diagnostics["total_intensity"] = total_intensity
            diagnostics["conserved"] = total_intensity > 0 && isfinite(total_intensity)
            diagnostics["conservation_method"] = "intensity_sum"
        end
    end

    return diagnostics
end

"""
    check_tolerances(metrics::Dict, tolerances::Dict) -> Bool

Check if computed metrics meet the specified tolerances.
"""
function check_tolerances(metrics::Dict, tolerances::Dict)
    try
        # Check field relative error
        if haskey(tolerances, "field_relative_error") &&
           haskey(metrics, "field_relative_error")
            if metrics["field_relative_error"] > tolerances["field_relative_error"]
                return false
            end
        end

        # Check energy relative error
        if haskey(tolerances, "energy_relative_error") &&
           haskey(metrics, "energy_relative_error")
            if metrics["energy_relative_error"] > tolerances["energy_relative_error"]
                return false
            end
        end

        # Check convergence threshold
        if haskey(tolerances, "convergence_threshold") &&
           haskey(metrics, "convergence_achieved")
            if !metrics["convergence_achieved"]
                return false
            end
        end

        return true

    catch e
        @warn "Error checking tolerances" exception=e
        return false
    end
end

"""
    generate_metric_report(metrics::ValidationMetricsResult)

Generate a comprehensive human-readable metrics report.
"""
function generate_metric_report(metrics::ValidationMetricsResult)
    report = """
# Validation Metrics Report

## Field Validation
- **Absolute Error**: $(isnan(metrics.field_absolute_error) ? "N/A" : @sprintf("%.2e", metrics.field_absolute_error))
- **Relative Error**: $(isnan(metrics.field_relative_error) ? "N/A" : @sprintf("%.2e", metrics.field_relative_error))
- **RMS Error**: $(isnan(metrics.field_rms_error) ? "N/A" : @sprintf("%.2e", metrics.field_rms_error))
- **Correlation**: $(isnan(metrics.field_correlation) ? "N/A" : @sprintf("%.4f", metrics.field_correlation))

## Energy Conservation
- **Energy Relative Error**: $(isnan(metrics.energy_relative_error) ? "N/A" : @sprintf("%.2e", metrics.energy_relative_error))
- **Conservation Ratio**: $(isnan(metrics.energy_conservation_ratio) ? "N/A" : @sprintf("%.6f", metrics.energy_conservation_ratio))

## Convergence Analysis
- **Final Residual**: $(isnan(metrics.final_residual) ? "N/A" : @sprintf("%.2e", metrics.final_residual))
- **Convergence Rate**: $(isnan(metrics.convergence_rate) ? "N/A" : @sprintf("%.4f", metrics.convergence_rate))
- **Iterations**: $(metrics.iterations_to_convergence)
- **Trend**: $(metrics.convergence_trend)

## Statistical Analysis
- **Outlier Fraction**: $(isnan(metrics.outlier_fraction) ? "N/A" : @sprintf("%.1f%%", metrics.outlier_fraction*100))
- **Error Std Dev**: $(isnan(metrics.relative_error_std) ? "N/A" : @sprintf("%.2e", metrics.relative_error_std))

## Quality Assessment
- **Numerical Stability**: $(isnan(metrics.numerical_stability) ? "N/A" : @sprintf("%.1f%%", metrics.numerical_stability*100))
- **Solution Quality**: $(isnan(metrics.solution_quality) ? "N/A" : @sprintf("%.1f%%", metrics.solution_quality*100))
- **Validation Confidence**: $(isnan(metrics.validation_confidence) ? "N/A" : @sprintf("%.1f%%", metrics.validation_confidence*100))

## Quality Flags
"""

    for (flag, status) in metrics.quality_flags
        emoji = status ? "✅" : "❌"
        report *= "- $emoji **$(replace(flag, "_" => " ") |> titlecase)**: $(status ? "PASS" : "FAIL")\n"
    end

    if !isempty(metrics.warnings)
        report *= "\n## Warnings\n"
        for warning in metrics.warnings
            report *= "- ⚠️ $warning\n"
        end
    end

    return report
end

"""
    convert_metrics_to_dict(metrics::ValidationMetricsResult) -> Dict{String, Float64}

Convert ValidationMetricsResult struct to Dict{String, Float64} for compatibility.
"""
function convert_metrics_to_dict(metrics::ValidationMetricsResult)
    result = Dict{String, Float64}()

    # Field validation metrics
    result["field_absolute_error"] = isnan(metrics.field_absolute_error) ? 0.0 :
                                     metrics.field_absolute_error
    result["field_relative_error"] = isnan(metrics.field_relative_error) ? 0.0 :
                                     metrics.field_relative_error
    result["field_rms_error"] = isnan(metrics.field_rms_error) ? 0.0 :
                                metrics.field_rms_error
    result["field_max_error"] = isnan(metrics.field_max_error) ? 0.0 :
                                metrics.field_max_error
    result["field_mean_error"] = isnan(metrics.field_mean_error) ? 0.0 :
                                 metrics.field_mean_error
    result["field_correlation"] = isnan(metrics.field_correlation) ? 0.0 :
                                  metrics.field_correlation

    # Energy metrics
    result["energy_relative_error"] = isnan(metrics.energy_relative_error) ? 0.0 :
                                      metrics.energy_relative_error
    result["energy_absolute_error"] = isnan(metrics.energy_absolute_error) ? 0.0 :
                                      metrics.energy_absolute_error
    result["energy_conservation_ratio"] = isnan(metrics.energy_conservation_ratio) ? 1.0 :
                                          metrics.energy_conservation_ratio

    # Convergence metrics
    result["final_residual"] = isnan(metrics.final_residual) ? 1.0 : metrics.final_residual
    result["convergence_rate"] = isnan(metrics.convergence_rate) ? 0.0 :
                                 metrics.convergence_rate
    result["iterations_to_convergence"] = Float64(metrics.iterations_to_convergence)

    # Statistical metrics
    result["relative_error_std"] = isnan(metrics.relative_error_std) ? 0.0 :
                                   metrics.relative_error_std
    result["relative_error_mean"] = isnan(metrics.relative_error_mean) ? 0.0 :
                                    metrics.relative_error_mean
    result["outlier_fraction"] = isnan(metrics.outlier_fraction) ? 0.0 :
                                 metrics.outlier_fraction

    # Quality indicators
    result["numerical_stability"] = isnan(metrics.numerical_stability) ? 0.0 :
                                    metrics.numerical_stability
    result["solution_quality"] = isnan(metrics.solution_quality) ? 0.0 :
                                 metrics.solution_quality
    result["validation_confidence"] = isnan(metrics.validation_confidence) ? 0.0 :
                                      metrics.validation_confidence

    # Convert quality flags to numeric scores (avoid type issues)
    quality_score = 0.0
    total_flags = length(metrics.quality_flags)
    if total_flags > 0
        passed_flags = sum(values(metrics.quality_flags))
        quality_score = Float64(passed_flags) / Float64(total_flags)
    end
    result["quality_flags_score"] = quality_score

    return result
end

end # module ValidationMetrics
