# Basic Scattering Example

This example demonstrates the fundamental workflow of FrequencyMaxwell.jl by simulating electromagnetic scattering from a spherical dielectric bead in water. This represents a common configuration in optical microscopy and particle analysis.

## Problem Setup

We simulate a 500 nm radius silica bead (n = 1.46) in water (n = 1.333) illuminated by a 532 nm plane wave. This configuration is representative of:

- Optical particle sizing
- Flow cytometry measurements
- Microscopic imaging of dielectric particles
- Light scattering spectroscopy

## Complete Example Code

```julia
using FrequencyMaxwell
using LinearSolve

function basic_scattering_example()
    println("FrequencyMaxwell Basic Scattering Example")
    println("=" ^ 45)

    # Step 1: Create electromagnetic solver with integrated configuration
    println("Creating solver with streamlined API...")
    solver = ConvergentBornSolver(
        wavelength = 532e-9,          # 532 nm (green laser)
        permittivity_bg = 1.333^2,    # Water background (n=1.333)
        resolution = (50e-9, 50e-9, 50e-9),  # 50 nm isotropic resolution
        grid_size = (128, 128, 64),   # 6.4 × 6.4 × 3.2 μm domain
        boundary_conditions = PeriodicBoundaryCondition(),  # Periodic boundaries
        tolerance = 1e-6,             # Convergence tolerance
        linear_solver = KrylovJL_GMRES()  # Recommended linear algebra algorithm
    )

    println("Configuration:")
    println("  Wavelength: $(solver.wavelength * 1e9) nm")
    println("  Background n: $(sqrt(solver.permittivity_bg))")
    println("  Grid size: $(solver.grid_size)")
    println("  Domain size: $(domain_size(solver) .* 1e6) μm")

    # Step 2: Define the incident plane wave source
    println("Setting up plane wave source...")
    source = PlaneWaveSource(
        wavelength = solver.wavelength,
        polarization = [1.0, 0.0, 0.0],  # X-polarized light
        k_vector = [0.0, 0.0, 1.0],      # Propagating in +Z direction
        amplitude = 1.0                   # 1 V/m amplitude
    )

    println("Source configuration:")
    println("  Wavelength: $(source_wavelength(source) * 1e9) nm")
    println("  Polarization: $(source.polarization)")
    println("  Propagation: $(source.k_vector)")

    # Step 3: Create a spherical bead phantom
    println("\\nGenerating phantom...")
    bead_radius_pixels = 10.0  # 500 nm radius (10 pixels × 50 nm)
    phantom = phantom_bead(
        solver.grid_size,
        [1.46^2],                 # Silica bead (n=1.46)
        bead_radius_pixels        # 500 nm radius
    )

    # Calculate phantom statistics
    bead_volume = count(real.(phantom) .> solver.permittivity_bg * 1.05)
    total_volume = prod(solver.grid_size)
    volume_fraction = bead_volume / total_volume

    println("Phantom statistics:")
    println("  Bead material: n = $(sqrt(1.46^2)) (silica)")
    println("  Bead radius: $(bead_radius_pixels * solver.resolution[1] * 1e9) nm")
    println("  Volume fraction: $(round(volume_fraction * 100, digits=2))%")

    # Step 4: Solve the electromagnetic scattering problem
    println("\\nSolving electromagnetic scattering...")
    EMfield = solve(solver, source, phantom)

    println("Solution completed successfully!")
    println("Field properties:")
    println("  Electric field size: $(size(EMfield.E))")
    println("  Magnetic field size: $(size(EMfield.H))")

    # Step 5: Analyze results
    println("\\nAnalyzing results...")

    # Calculate field intensity
    intensity = field_intensity(EMfield)
    max_intensity = maximum(intensity)
    mean_intensity = sum(intensity) / length(intensity)
    enhancement = max_intensity / 1.0  # Relative to incident intensity

    println("  Maximum intensity: $(round(max_intensity, digits=4)) (V/m)²")
    println("  Mean intensity: $(round(mean_intensity, digits=4)) (V/m)²")
    println("  Enhancement factor: $(round(enhancement, digits=2))")

    # Calculate total field energy
    total_energy = field_energy(EMfield)

    println("  Total field energy: $(total_energy) J")

    # Extract central plane for analysis
    z_center = div(solver.grid_size[3], 2)
    central_plane = extract_plane(EMfield, 3, z_center)
    plane_intensity = field_intensity(central_plane)

    println("Central plane analysis:")
    println("  Plane size: $(size(plane_intensity))")
    println("  Plane mean intensity: $(round(sum(plane_intensity)/length(plane_intensity), digits=4))")

    println("\\nExample completed successfully!")
    return EMfield, phantom
end

# Run the example
EMfield, phantom = basic_scattering_example()
```

## Step-by-Step Explanation

### 1. Solver Configuration

The `ConvergentBornSolver` defines all physical and numerical parameters:

```julia
solver = ConvergentBornSolver(
    wavelength = 532e-9,          # Physical wavelength in vacuum
    permittivity_bg = 1.333^2,    # Background relative permittivity
    resolution = (50e-9, 50e-9, 50e-9),  # Voxel size
    grid_size = (128, 128, 64),   # Number of grid points
    boundary_conditions = PeriodicBoundaryCondition(),  # Boundary conditions
    tolerance = 1e-6,             # Convergence criterion
    linear_solver = KrylovJL_GMRES()  # Recommended linear algebra algorithm
)
```

**Key considerations**:
- **Wavelength**: Always specify vacuum wavelength (532 nm for green laser)
- **Resolution**: 50 nm provides ~10 points per wavelength in water
- **Grid size**: Large enough to contain the object with adequate boundary padding
- **Background permittivity**: εᵣ = n² = 1.333² ≈ 1.77 for water
- **Boundary conditions**: Use `PeriodicBoundaryCondition()` for periodic structures or `AbsorbingBoundaryCondition()` for open domains

### 2. Electromagnetic Solving

The solver object contains:
- All configuration parameters directly accessible (e.g., `solver.wavelength`)
- Preallocated arrays for efficient field computation
- Green's function operators
- LinearSolve.jl integration components

### 3. Source Definition

```julia
source = PlaneWaveSource(
    wavelength = solver.wavelength,
    polarization = [1.0, 0.0, 0.0],  # Linear polarization
    k_vector = [0.0, 0.0, 1.0],      # Propagation direction
    amplitude = 1.0                   # Field amplitude
)
```

**Important points**:
- Polarization must be perpendicular to k_vector (automatically validated)
- k_vector is automatically normalized to unit magnitude
- Amplitude sets the incident field strength

### 4. Phantom Generation

```julia
phantom = phantom_bead(
    solver.grid_size,
    [1.46^2],         # Relative permittivity of silica
    10.0              # Radius in grid points
)
```

This creates a 3D array with:
- Background regions: permittivity = 1.0 (will be scaled by solver.permittivity_bg)
- Bead region: permittivity = 1.46² ≈ 2.13

### 5. Electromagnetic Solving

```julia
EMfield = solve(solver, source, phantom)
```

The solver:
1. Generates incident fields from the source
2. Sets up the linear system (I - G·V)E = G·E_inc
3. Solves iteratively using the specified linear algebra algorithm
4. Computes total fields (incident + scattered)

### 6. Results Analysis

The `ElectromagneticField` structure contains:

```julia
# Field arrays with dimensions (nx, ny, nz, 3)
EMfield.E  # Electric field vector components
EMfield.H  # Magnetic field vector components

# Analysis functions
intensity = field_intensity(EMfield)      # |E|² intensity
energy = field_energy(EMfield)            # ∫ |E|² dV
poynting = poynting_vector(EMfield)       # Power flow
```

## Physical Interpretation

### Scattering Physics

The silica bead (n = 1.46) in water (n = 1.333) creates a refractive index contrast that scatters the incident light. Key physics:

- **Forward scattering**: Concentrated in the original propagation direction
- **Field enhancement**: Increased intensity near the particle surface
- **Phase shifts**: Complex field patterns due to optical path differences

### Parameter Scaling

This example uses realistic optical parameters:

- **Particle size**: 500 nm radius ≈ wavelength in water (λ/n ≈ 400 nm)
- **Size parameter**: ka ≈ 7.8 (intermediate scattering regime)
- **Refractive contrast**: Δn ≈ 0.13 (moderate contrast)

### Expected Results

For this configuration, expect:
- **Enhancement factor**: 2-4× intensity enhancement near the particle
- **Convergence**: 20-50 iterations for 10⁻⁶ tolerance
- **Computation time**: Seconds to minutes depending on hardware

## Variations and Extensions

### Different Materials

```julia
# High contrast: polystyrene bead
phantom_ps = phantom_bead(solver.grid_size, [1.59^2], 10.0)

# Metallic particle (with loss)
n_gold = 0.47 + 2.4im  # Gold at 532 nm
phantom_gold = phantom_bead(solver.grid_size, [n_gold^2], 5.0)

# Biological cell
phantom_cell = phantom_bead(solver.grid_size, [1.38^2], 15.0)
```

### Different Illumination

```julia
# Oblique incidence
source_oblique = PlaneWaveSource(
    wavelength = 532e-9,
    polarization = [0.0, 1.0, 0.0],  # s-polarized
    k_vector = [sin(30°), 0.0, cos(30°)],  # 30° incidence
    amplitude = 1.0
)

# Circular polarization
source_circular = PlaneWaveSource(
    wavelength = 532e-9,
    polarization = [1.0, 1im, 0.0],  # Right circular
    k_vector = [0.0, 0.0, 1.0],
    amplitude = 1.0
)
```

### Multiple Particles

```julia
# Create two beads
phantom_dual = ones(ComplexF64, solver.grid_size)

# First bead at offset position
phantom1 = phantom_bead(solver.grid_size, [1.46^2], 8.0)
# Apply offset by manipulating indices...

# Second bead at different position
phantom2 = phantom_bead(solver.grid_size, [1.59^2], 6.0)
# Apply different offset...

# Combine phantoms
phantom_dual = phantom1 .* (phantom2 .== 1.0) + phantom2 .* (phantom2 .!= 1.0)
```

### Performance Optimization

```julia
# Use Float32 for memory efficiency on large problems
solver_optimized = ConvergentBornSolver(
    wavelength = 532e-9,
    permittivity_bg = 1.333^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (256, 256, 128),  # Larger problem
    boundary_conditions = PeriodicBoundaryCondition(),
    tolerance = 1e-6,
    linear_solver = KrylovJL_GMRES()  # Recommended default
)
# Optimized for large-scale problems
```

## Troubleshooting

### Convergence Issues

If the solver doesn't converge:

```julia
# Try different linear solver
solver = ConvergentBornSolver(
    wavelength = 532e-9,
    permittivity_bg = 1.333^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (128, 128, 64),
    boundary_conditions = PeriodicBoundaryCondition(),
    linear_solver = KrylovJL_BICGSTAB(),  # Alternative algorithm
    tolerance = 1e-5               # Relax tolerance
)
```

### Memory Issues

For large problems:

```julia
# Reduce grid size during development
solver = ConvergentBornSolver(
    wavelength = 532e-9,
    permittivity_bg = 1.333^2,
    resolution = (50e-9, 50e-9, 50e-9),
    grid_size = (64, 64, 32),  # Smaller for testing
    boundary_conditions = PeriodicBoundaryCondition(),
    tolerance = 1e-6
)
```

### Physical Validation

Check results for physical consistency:

```julia
# Verify energy conservation
incident_power = source_power(source)
scattered_power = field_energy(EMfield) - incident_power
absorption = # ... calculate if materials are lossy

println("Power balance check:")
println("  Incident: $incident_power")
println("  Scattered: $scattered_power")
```

This basic scattering example provides the foundation for more complex electromagnetic simulations using FrequencyMaxwell.jl.