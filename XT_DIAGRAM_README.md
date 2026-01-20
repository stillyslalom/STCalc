# X-T Diagram Calculator

A 1D multicomponent compressible Euler solver for generating shock tube x-t diagrams using JavaScript running in your web browser.

## Overview

The X-T Diagram Calculator solves the 1D compressible Euler equations using a finite volume method with HLLC (Harten-Lax-van Leer-Contact) flux to simulate shock tube evolution. It generates space-time (x-t) diagrams showing pressure fields and interface trajectories.

## Features

### Core Capabilities
- **Multicomponent Gas Support**: Configure shock tubes with multiple gas slabs of different species
- **HLLC Flux Scheme**: High-resolution Riemann solver for accurate shock capturing
- **Lagrangian Tracers**: Track material interfaces through space and time
- **Interactive Visualization**: Hover over the x-t diagram to see pressure values
- **Data Export**: Download results as CSV or JSON for further analysis

### Gas Database
- Pre-defined gases: Air, Nitrogen, Helium, Argon, Xenon, CO₂, SF₆, Neon, Acetone Vapor
- Custom gas mixtures supported
- Import initial conditions from the main Shock Tube Calculator

## Technical Details

### Numerical Method
- **Spatial Discretization**: Finite volume method on uniform grid
- **Flux Function**: HLLC Riemann solver with acoustic wave speed estimates
- **Time Integration**: 2nd-order Runge-Kutta (RK2) method
- **CFL Condition**: Adaptive time stepping with configurable CFL number (default 0.4)
- **Boundary Conditions**: Reflective walls at domain boundaries

### Governing Equations
The solver implements the 1D Euler equations in conservative form:

```
∂U/∂t + ∂F/∂x = 0
```

Where:
- U = [ρ, ρu, E]ᵀ (density, momentum, total energy)
- F = [ρu, ρu² + p, u(E + p)]ᵀ (flux vector)
- p = (γ-1)(E - ½ρu²) (equation of state for ideal gas)

### Gas Properties
Each gas region is characterized by:
- **γ (gamma)**: Specific heat ratio
- **MW**: Molecular weight (g/mol)
- **R**: Gas constant = Rᵤ/MW where Rᵤ = 8314.51 J/(kmol·K)

### Performance
- Typical grid: 500 spatial points
- Simulation time: ~20 ms of physical time
- Computation time: 1-3 seconds for standard cases
- X-T snapshots: ~200 time levels stored for visualization

## Usage

### Basic Workflow

1. **Configure Gas Slabs**
   - Click "Add Gas Slab" to add new sections
   - Specify: Gas type, Pressure (kPa), Temperature (K), Length (m)
   - Arrange slabs in order using ↑↓ buttons
   - Remove unwanted slabs with "Remove" button

2. **Import from Calculator** (Optional)
   - Run a calculation in the main Shock Tube Calculator
   - Click "Import from Shock Tube Calculator"
   - Adjust slab lengths as needed

3. **Set Simulation Parameters**
   - **Grid Points**: 100-2000 (more = higher resolution, slower computation)
   - **Simulation Time**: Duration in milliseconds
   - **CFL Number**: Stability parameter (0.1-0.9, recommended: 0.4)

4. **Run Simulation**
   - Click "Run Simulation"
   - Watch progress bar
   - View x-t diagram and statistics

5. **Export Results**
   - **Download Image (PNG)**: Export the visualization
   - **Download X-T Data (CSV)**: Full pressure field at all time/space points
   - **Download Tracer Data (CSV)**: Interface trajectory coordinates
   - **Download Complete Data (JSON)**: All data including configuration

### Example Configuration

**Classic Helium-Air Shock Tube:**
- Slab 1: Air, 101 kPa, 300 K, 6 m (driven section)
- Slab 2: Helium, 400 kPa, 300 K, 3 m (driver section)
- Grid Points: 500
- Simulation Time: 20 ms
- CFL: 0.4

## Data Export Formats

### X-T Data CSV
```csv
Time (ms),Position (m),Pressure (kPa)
0.0000,0.0009,101.33
0.0000,0.0027,101.33
...
```

### Tracer Data CSV
```csv
Interface,Time (ms),Position (m)
1,0.0000,6.0000
1,0.1000,6.0531
...
```

### Complete JSON
```json
{
  "configuration": {
    "slabs": [...],
    "gridPoints": 500,
    "simulationTime": 20,
    "cfl": 0.4
  },
  "results": {
    "x": [...],
    "t": 0.02,
    "p": [...],
    "u": [...],
    "rho": [...],
    "T": [...],
    "xtData": [...],
    "tracers": [...]
  }
}
```

## Visualization Features

### X-T Diagram
- **Heatmap**: Pressure field colored using perceptually-uniform plasma colormap
- **White Lines**: Interface trajectories showing material discontinuities
- **Axes**: Position (m) on x-axis, Time (ms) on y-axis
- **Colorbar**: Pressure scale in kPa

### Interactive Features
- **Hover**: Display exact pressure, position, and time values
- **Crosshairs**: Visual guides when hovering

## Algorithm Details

### HLLC Flux Calculation
The HLLC solver resolves the Riemann problem at each cell interface:

1. Compute left and right states (UL, UR)
2. Calculate wave speeds: SL, SR, S* (left, right, and contact)
3. Determine flux based on wave structure:
   - If SL ≥ 0: Use left flux
   - If SL < 0 < S*: Use left star flux
   - If S* < 0 < SR: Use right star flux
   - If SR ≤ 0: Use right flux

### Lagrangian Tracers
Interface positions are tracked by integrating:
```
dx/dt = u(x,t)
```
using forward Euler with the local fluid velocity.

## Limitations

1. **1D Assumption**: No transverse effects or boundary layers
2. **Ideal Gas**: Assumes calorically perfect gas (constant γ)
3. **Inviscid**: No viscosity or heat conduction
4. **No Chemistry**: No reactions or dissociation
5. **Reflective Boundaries**: Tube ends treated as solid walls

## Browser Compatibility

Tested on:
- Chrome 90+
- Firefox 88+

Requires:
- HTML5 Canvas support
- JavaScript ES6+
- LocalStorage for importing calculator settings

## Performance Tips

1. **Start Small**: Begin with 200-300 grid points for quick testing
2. **Increase Resolution**: Use 500-1000 points for publication-quality results
3. **Simulation Time**: Match physical phenomena (shocks reach walls ~10-20 ms typically)
4. **CFL Number**: Lower values (0.3-0.4) are more stable but slower

## References

1. Toro, E.F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics*. Springer.
2. LeVeque, R.J. (2002). *Finite Volume Methods for Hyperbolic Problems*. Cambridge University Press.
3. Anderson, J.D. (2003). *Modern Compressible Flow*. McGraw-Hill.

## Implementation Notes

### File Structure
```
js/
├── euler-solver.js       # Core solver with HLLC flux
├── xt-visualization.js   # Canvas-based visualization
└── shock-import.js       # Calculator integration

css/
└── xt-diagram.css        # Styling

xt-diagram.html           # Main application page
```

### Key Classes

**EulerSolver**
- Properties: nx, cfl, finalTime, U, flux, tracers, xtData
- Methods: initialize(), step(), run(), computeFluxes(), hllcFlux()

**XTVisualization**
- Properties: canvas, ctx, xtData, tracers, colormap
- Methods: render(), drawHeatmap(), drawTracers(), plasmaColormap()

## Contributing

Issues and feature requests: https://github.com/stillyslalom/STCalc

## License

Same as parent project (check repository for details).

## Authors

Alex Ames, 2026

Based on the original Shock Tube Calculator by Jason Oakley.
