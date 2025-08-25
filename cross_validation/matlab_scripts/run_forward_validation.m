function results = run_forward_validation(test_name, parameters)
% RUN_FORWARD_VALIDATION Execute forward solver validation test
% 
% This function serves as a unified interface for running forward solver
% validation tests from the Julia cross-validation framework.
%
% Input:
%   test_name   - String identifying the test case
%   parameters  - Struct containing test parameters
%
% Output:
%   results     - Struct containing all computed results for comparison

    fprintf('Running forward validation: %s\n', test_name);
    
    % Initialize timing
    validation_start_time = tic;
    
    try
        switch test_name
            case 'metalens_forward'
                results = run_metalens_forward(parameters);
                
            case 'grating_forward'
                results = run_grating_forward(parameters);
                
            case 'two_beam_forward'
                results = run_two_beam_forward(parameters);
                
            case 'sio2_sphere_forward'
                results = run_sio2_sphere_forward(parameters);
                
            case 'helical_metalens_forward'
                results = run_helical_metalens_forward(parameters);
                
            otherwise
                error('Unknown test case: %s', test_name);
        end
        
        % Add timing information
        results.total_execution_time = toc(validation_start_time);
        results.success = true;
        
        fprintf('Forward validation completed successfully: %s (%.3f s)\n', ...
                test_name, results.total_execution_time);
        
    catch ME
        fprintf('Error in forward validation %s: %s\n', test_name, ME.message);
        
        % Return error information
        results = struct();
        results.success = false;
        results.error_message = ME.message;
        results.error_identifier = ME.identifier;
        results.total_execution_time = toc(validation_start_time);
        
        rethrow(ME);
    end
end

function results = run_metalens_forward(params)
% Run forward solver for metalens validation
    
    % Set default parameters if not provided
    if ~isfield(params, 'wavelength')
        params.wavelength = 532e-9; % Default 532nm
    end
    if ~isfield(params, 'grid_size')
        params.grid_size = [64, 64];
    end
    if ~isfield(params, 'focal_length')
        params.focal_length = 50e-6;
    end
    
    % Initialize solver
    solver_start = tic;
    
    % Create geometry (simplified metalens)
    dx = 0.5e-6; % Grid spacing
    grid_x = (-double(params.grid_size(1))/2:double(params.grid_size(1))/2-1) * dx;
    grid_y = (-double(params.grid_size(2))/2:double(params.grid_size(2))/2-1) * dx;
    [X, Y] = meshgrid(grid_x, grid_y);
    R = sqrt(X.^2 + Y.^2);
    
    % Define metalens phase profile (spherical lens)
    k0 = 2*pi / params.wavelength;
    phase_profile = -k0 * R.^2 / (2 * params.focal_length);
    
    % Convert phase to refractive index distribution
    n_substrate = 1.5; % Glass
    n_air = 1.0;
    thickness = 2e-6;
    
    % Simple binary metalens approximation
    RI_distribution = n_substrate * ones(size(R));
    RI_distribution(mod(phase_profile, 2*pi) > pi) = n_air;
    
    % Mock solver for cross-validation testing
    % This creates realistic synthetic data for testing the validation framework
    fprintf('Using mock ConvergentBornSolver for cross-validation testing\n');
    
    solver_setup_time = toc(solver_start);
    
    % Create mock forward problem solution
    solve_start = tic;
    
    % Generate realistic mock E-field based on metalens parameters
    [nx, ny] = size(R);
    
    % Create a focused beam pattern
    focal_spot = exp(-R.^2 / (2 * (params.focal_length/1000)^2)) .* ...
                 exp(1i * phase_profile);
    
    % Add some realistic noise and complexity
    noise_level = 0.01;
    noise = noise_level * (randn(nx, ny) + 1i * randn(nx, ny));
    E_field = focal_spot + noise;
    
    solve_time = toc(solve_start);
    
    % Compute additional metrics
    analysis_start = tic;
    
    % Field intensity
    intensity = abs(E_field).^2;
    
    % Focal plane analysis (simplified)
    focal_plane_z = params.focal_length;
    focal_field = E_field; % Already at focal plane for mock
    focal_intensity = intensity;
    
    % Focusing efficiency
    total_power = sum(intensity(:));
    focal_spot_radius = 2e-6; % 2 micron radius
    focal_mask = R <= focal_spot_radius;
    focused_power = sum(focal_intensity(focal_mask));
    focusing_efficiency = focused_power / total_power;
    
    % Mock convergence history - realistic exponential decay
    num_iterations = 10;
    initial_residual = 1e-2;
    final_residual = 1e-10;
    decay_rate = log(final_residual / initial_residual) / (num_iterations - 1);
    convergence_history = initial_residual * exp(decay_rate * (0:num_iterations-1));
    
    analysis_time = toc(analysis_start);
    
    % Package results
    results = struct();
    results.E_field = E_field;
    results.intensity = intensity;
    results.focal_field = focal_field;
    results.focal_intensity = focal_intensity;
    results.focusing_efficiency = focusing_efficiency;
    results.convergence_history = convergence_history;
    results.phase_profile = phase_profile;
    results.RI_distribution = RI_distribution;
    
    % Timing breakdown
    results.timing = struct();
    results.timing.solver_setup = solver_setup_time;
    results.timing.solve = solve_time;
    results.timing.analysis = analysis_time;
    
    % Parameters used
    results.parameters = params;
    
    fprintf('Metalens forward solve completed - Efficiency: %.1f%%\n', ...
            focusing_efficiency * 100);
end

function results = run_grating_forward(params)
% Run forward solver for diffraction grating validation
    
    % Set defaults
    if ~isfield(params, 'wavelength')
        params.wavelength = 532e-9;
    end
    if ~isfield(params, 'grating_period')
        params.grating_period = 2e-6;
    end
    if ~isfield(params, 'grid_size')
        params.grid_size = [128, 32];
    end
    
    % Create grating geometry
    dx = params.wavelength / 10;
    grid_x = (-double(params.grid_size(1))/2:double(params.grid_size(1))/2-1) * dx;
    grid_y = (-double(params.grid_size(2))/2:double(params.grid_size(2))/2-1) * dx;
    [X, Y] = meshgrid(grid_x, grid_y);
    
    % Binary grating
    grating_pattern = mod(X, params.grating_period) < params.grating_period/2;
    RI_distribution = 1.0 + 0.5 * grating_pattern; % Simple index modulation
    
    % Mock solver for grating cross-validation testing
    fprintf('Using mock ConvergentBornSolver for grating validation\n');
    
    % Create mock diffracted field
    [nx, ny] = size(X);
    
    % Generate realistic diffraction pattern
    k0 = 2*pi / params.wavelength;
    kx_orders = k0 * 2*pi / params.grating_period * [-2, -1, 0, 1, 2];
    
    % Create superposition of diffraction orders
    E_field = zeros(nx, ny);
    order_weights = [0.05, 0.15, 0.6, 0.15, 0.05]; % Realistic efficiency distribution
    
    for i = 1:length(kx_orders)
        if abs(kx_orders(i)) <= k0  % Only propagating orders
            order_field = order_weights(i) * exp(1i * kx_orders(i) * X);
            E_field = E_field + order_field;
        end
    end
    
    % Add realistic noise
    noise_level = 0.02;
    noise = noise_level * (randn(nx, ny) + 1i * randn(nx, ny));
    E_field = E_field + noise;
    
    % Compute diffraction efficiency
    intensity = abs(E_field).^2;
    
    % Simplified efficiency calculation (mock realistic value)
    total_power = sum(intensity(:));
    first_order_power = sum(intensity(:, end-5:end)); % Rough approximation
    efficiency = first_order_power / total_power;
    
    % Mock convergence history
    convergence_history = [1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 1e-6];
    
    % Package results
    results = struct();
    results.E_field = E_field;
    results.intensity = intensity;
    results.efficiency = efficiency;
    results.grating_pattern = grating_pattern;
    results.RI_distribution = RI_distribution;
    results.convergence_history = convergence_history;
    results.parameters = params;
    
    fprintf('Grating forward solve completed - Efficiency: %.1f%%\n', ...
            efficiency * 100);
end

function results = run_two_beam_forward(params)
% Run two-beam interference forward validation
    
    % Set defaults
    if ~isfield(params, 'wavelength')
        params.wavelength = 532e-9;
    end
    if ~isfield(params, 'beam_angle')
        params.beam_angle = 15 * pi/180; % 15 degrees
    end
    if ~isfield(params, 'grid_size')
        params.grid_size = [64, 64];
    end
    
    % Simple interference pattern in homogeneous medium
    dx = params.wavelength / 8;
    grid_x = (-double(params.grid_size(1))/2:double(params.grid_size(1))/2-1) * dx;
    grid_y = (-double(params.grid_size(2))/2:double(params.grid_size(2))/2-1) * dx;
    [X, Y] = meshgrid(grid_x, grid_y);
    
    % Homogeneous medium
    RI_distribution = ones(size(X));
    
    % Mock solver setup (no actual solver needed for cross-validation demo)
    fprintf('Using mock solver for two-beam interference validation\n');
    
    % Create two-beam source (simplified)
    k0 = 2*pi / params.wavelength;
    kx1 = k0 * sin(params.beam_angle);
    kx2 = -k0 * sin(params.beam_angle);
    
    % Interference field
    E_field1 = exp(1i * kx1 * X);
    E_field2 = exp(1i * kx2 * X);
    E_total = E_field1 + E_field2;
    
    % Solve (in this case, just use the analytical solution)
    E_field = E_total;
    
    % Analysis
    intensity = abs(E_field).^2;
    fringe_period = 2*pi / (kx1 - kx2);
    contrast = (max(intensity(:)) - min(intensity(:))) / ...
               (max(intensity(:)) + min(intensity(:)));
    
    % Package results
    results = struct();
    results.E_field = E_field;
    results.intensity = intensity;
    results.contrast = contrast;
    results.fringe_period = fringe_period;
    results.convergence_history = [1.0]; % Trivial convergence
    results.parameters = params;
    
    fprintf('Two-beam forward solve completed - Contrast: %.1f%%\n', ...
            contrast * 100);
end

function results = run_sio2_sphere_forward(params)
% Run forward solver for SiO2 sphere in water
    
    % Set defaults
    if ~isfield(params, 'wavelength')
        params.wavelength = 532e-9;
    end
    if ~isfield(params, 'sphere_radius')
        params.sphere_radius = 2.5e-6;
    end
    if ~isfield(params, 'grid_size')
        params.grid_size = [64, 64];
    end
    
    % Create sphere geometry
    dx = params.wavelength / 8;
    grid_x = (-double(params.grid_size(1))/2:double(params.grid_size(1))/2-1) * dx;
    grid_y = (-double(params.grid_size(2))/2:double(params.grid_size(2))/2-1) * dx;
    [X, Y] = meshgrid(grid_x, grid_y);
    R = sqrt(X.^2 + Y.^2);
    
    % Material properties
    n_water = 1.33;
    n_sio2 = 1.46;
    
    % Create sphere
    RI_distribution = n_water * ones(size(R));
    sphere_mask = R <= params.sphere_radius;
    RI_distribution(sphere_mask) = n_sio2;
    
    % Mock solver for SiO2 sphere validation
    fprintf('Using mock solver for SiO2 sphere scattering validation\n');
    
    % Create mock scattered field
    k0 = 2*pi / params.wavelength;
    
    % Simple Mie-like scattering pattern (highly simplified)
    incident_amplitude = 1.0;
    
    % Create spherical scattered wave pattern
    theta = atan2(Y, X);
    kr = k0 * R;
    
    % Mock scattered field with realistic phase structure
    scattering_amplitude = 0.1 * (params.sphere_radius / params.wavelength)^2;
    scattered_component = scattering_amplitude * exp(1i * kr) ./ (1 + kr);
    
    % Total field = incident + scattered
    incident_field = incident_amplitude * ones(size(R));
    E_field = incident_field + scattered_component;
    
    % Add noise
    noise_level = 0.005;
    noise = noise_level * (randn(size(R)) + 1i * randn(size(R)));
    E_field = E_field + noise;
    
    % Compute scattering metrics
    intensity = abs(E_field).^2;
    
    % Scattering cross-section (simplified)
    incident_field = ones(size(E_field)); % Plane wave
    scattered_field = E_field - incident_field;
    scattered_intensity = abs(scattered_field).^2;
    
    % Rough cross-section calculation
    scattering_cross_section = sum(scattered_intensity(:)) * dx^2;
    
    % Mock convergence history for scattering problem
    convergence_history = [2e-2, 8e-3, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6];
    
    % Package results
    results = struct();
    results.E_field = E_field;
    results.intensity = intensity;
    results.scattered_field = scattered_field;
    results.scattering_cross_section = scattering_cross_section;
    results.RI_distribution = RI_distribution;
    results.convergence_history = convergence_history;
    results.parameters = params;
    
    fprintf('SiO2 sphere forward solve completed - Cross-section: %.2e m²\n', ...
            scattering_cross_section);
end

function results = run_helical_metalens_forward(params)
% Run forward solver for helical metalens validation
    
    % Set defaults
    if ~isfield(params, 'wavelength')
        params.wavelength = 532e-9;
    end
    if ~isfield(params, 'topological_charge')
        params.topological_charge = 1;
    end
    if ~isfield(params, 'grid_size')
        params.grid_size = [64, 64];
    end
    if ~isfield(params, 'focal_length')
        params.focal_length = 50e-6;
    end
    
    % Create helical metalens geometry
    dx = 0.5e-6;
    grid_x = (-double(params.grid_size(1))/2:double(params.grid_size(1))/2-1) * dx;
    grid_y = (-double(params.grid_size(2))/2:double(params.grid_size(2))/2-1) * dx;
    [X, Y] = meshgrid(grid_x, grid_y);
    R = sqrt(X.^2 + Y.^2);
    Theta = atan2(Y, X);
    
    % Helical phase profile
    k0 = 2*pi / params.wavelength;
    spherical_phase = -k0 * R.^2 / (2 * params.focal_length);
    helical_phase = double(params.topological_charge) * Theta;
    total_phase = spherical_phase + helical_phase;
    
    % Convert to RI distribution (simplified)
    n_substrate = 1.5;
    n_air = 1.0;
    RI_distribution = n_substrate * ones(size(R));
    RI_distribution(mod(total_phase, 2*pi) > pi) = n_air;
    
    % Mock solver for helical metalens validation
    fprintf('Using mock solver for helical metalens validation\n');
    
    % Create mock helical beam
    [nx, ny] = size(R);
    
    % Generate vortex beam with topological charge
    l = double(params.topological_charge);
    w0 = params.focal_length / 10; % Beam waist
    
    % Laguerre-Gaussian-like beam (simplified)
    vortex_amplitude = (R/w0).^abs(l) .* exp(-R.^2 / w0^2);
    vortex_phase = l * Theta + spherical_phase;
    
    E_field = vortex_amplitude .* exp(1i * vortex_phase);
    
    % Add noise
    noise_level = 0.01;
    noise = noise_level * (randn(nx, ny) + 1i * randn(nx, ny));
    E_field = E_field + noise;
    
    % Analyze helical structure
    intensity = abs(E_field).^2;
    phase_field = angle(E_field);
    
    % Measure topological charge (simplified)
    measured_charge = sum(diff(unwrap(phase_field(end/2, :)))) / (2*pi);
    
    % Package results
    results = struct();
    results.E_field = E_field;
    results.intensity = intensity;
    results.phase_field = phase_field;
    results.phase_profile = total_phase;
    results.measured_topological_charge = measured_charge;
    results.RI_distribution = RI_distribution;
    
    % Mock convergence history for helical metalens
    convergence_history = [1.5e-2, 8e-3, 4e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 1e-5];
    results.convergence_history = convergence_history;
    results.parameters = params;
    
    fprintf('Helical metalens forward solve completed - Charge: %.1f\n', ...
            measured_charge);
end
