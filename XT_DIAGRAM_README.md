# x-t diagram

A 1D multicomponent compressible Euler solver that runs in the browser and draws shock tube x-t diagrams.

## Using it

1. Open **x-t diagram** on the page. If the shock tube calculation has been run, the driver and driven slabs start from its states; otherwise a helium/air example is shown. **Use driver and driven states from the shock tube calculation** re-imports them at any time and leaves slabs 3 and up alone.
2. Edit the slabs. They are listed from the driver end wall to the driven end wall. Each has a gas, liquid or mixture, a pressure, a temperature and a length. Add more slabs for test-section gases, liquids or buffer layers.
3. Optionally enter probe positions to get time histories at fixed locations, as a pressure transducer would record.
4. Run. Choose the plotted field (pressure, density, temperature, velocity) and a linear or logarithmic color scale afterwards without re-running. Hover over the diagram to read all four values at a point. Pressure and temperature use the units chosen at the top of the page; exported files always use kPa and K.
5. Export the image, the full x-t field, interface paths or probe traces as CSV, or everything as JSON or HDF5.

## Numerical method

- Finite volume on a uniform grid with the HLLC flux and Davis wave-speed estimates.
- Second order in space by default (MUSCL): density, velocity, pressure and volume fractions are reconstructed linearly in each cell with the minmod limiter, and the energy is rebuilt from them at each face. First order is available as an option.
- Time integration: SSP Runge–Kutta SSPRK(4,3) (default, CFL 0.8) or two-stage Runge–Kutta (CFL 0.4). Time steps are shortened to land exactly on snapshot times.
- Reflecting walls at both ends.
- Lagrangian tracers follow each initial interface and are drawn as white lines.

### Equation of state

Each slab has its own equation of state, and all materials in a cell share one pressure.

- **Constant γ** (default): a calorically perfect gas, or a stiffened gas for liquids:

  ```
  p = (γ − 1) ρ e − γ p∞        p + p∞ = ρ R T        a² = γ (p + p∞) / ρ
  ```

  `p∞ = 0` gives an ideal gas. The built-in water uses γ = 6.12, p∞ = 343 MPa, with R chosen to give 998.2 kg/m³ at 293.15 K and 1 atm (sound speed 1450 m/s). The temperature of a liquid is only as meaningful as that fit.

- **Variable γ**: gases are ideal, `p = ρ R T`, with specific heats from NASA polynomial fits and frozen composition. The internal energy `e(T)` of each slab's gas is tabulated every 2 K. Liquids keep the stiffened-gas model. About 1.5 times slower.

### Interfaces

- **Mixed cell** (default): a volume fraction per slab is advected with the flow using the HLLC contact velocity, inside the Runge–Kutta stages, which keeps pressure and velocity uniform across a moving material interface. With constant γ, `1/(γ−1)` and `γp∞/(γ−1)` are volume-averaged in mixed cells. With variable γ, each slab's partial density is also carried, so every material has its own density and temperature in a mixed cell, and the shared pressure is found by Newton iteration. Required for gas–liquid interfaces and for variable γ.
- **Sharp**: each cell takes the properties of the region between tracers. Faster, but not conservative at interfaces (about 1% error in wave strengths) and limited to slabs with the same p∞.

### Limits

- Frozen composition: no dissociation, ionization, reaction or phase change (no cavitation). With constant γ there is no vibrational excitation either.
- Inviscid and one-dimensional: no boundary layers, diaphragm opening or area changes.
- Shocks are spread over two or three cells and contacts over more (about 7 cells after 200 steps at second order, twice that at first order); increase the grid cells to sharpen them.

## Output files

- **x-t data (CSV)**: one row per cell and snapshot: time, position, density, velocity, pressure, temperature, γ, molecular weight, p∞, material.
- **Interface paths (CSV)**: position and velocity of each interface at every time step.
- **Probe traces (CSV)**: pressure, density, temperature and velocity at each probe at every snapshot.
- **HDF5**: groups `metadata`, `species` (γ, molecular weight, stiffening pressure, names), `coordinates` (`x`, `t`), `fields` (`density`, `velocity`, `pressure`, `temperature`, `gamma`, `molecular_weight`, `volume_fraction`, shaped snapshots × cells) and `tracers`. The HDF5 library is downloaded the first time this export is used.
- **JSON**: the run configuration and all results.

## Validation

`npm test` checks the solver against exact Riemann solutions (helium/air shock tube, air shock reflecting from water), checks the variable-γ solver against the variable-γ shock tube calculation, and checks that a moving air/water interface stays at uniform pressure and velocity with either gas model.
