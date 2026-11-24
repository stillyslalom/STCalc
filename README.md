This is the venerable [WiSTL shock tube calculator](https://silver.neep.wisc.edu/~shock/tools/gdcalc.html) with some updates. Credit is due to Jason Oakley & Riccardo Bonazza for the original calculator, which has proved very useful for shock tube work over the years.

Hosted at [https://stillyslalom.github.io/STCalc](https://stillyslalom.github.io/STCalc)

![Shock tube calculator screenshot](input_screenshot.png)

## What it does

- **Shock tube**: incident and reflected shock states and the driver pressure for a desired Mach number, or the Mach number produced by a given driver pressure. Unreachable Mach numbers and results outside the constant-γ range produce warnings.
- **Custom gases, liquids and mixtures**: define your own gases, mix gases by mole fraction, and model liquids such as water with the stiffened-gas equation of state. Available in every gas list.
- **Shocked interface (RMI)**: exact Riemann solution for the shock crossing an interface into a second gas or liquid, with Richtmyer–Meshkov growth estimates.
- **x-t diagram**: in-browser 1D multicomponent Euler solver with pressure, density, temperature and velocity maps, probe traces, and CSV/JSON/HDF5 export. See [XT_DIAGRAM_README.md](XT_DIAGRAM_README.md).
- **Timing** and **boost tank** helpers.
- **Shareable links**: "Copy link to this setup" encodes the shock tube inputs, including custom gases, in the URL.

## Layout

| Path | Contents |
| --- | --- |
| `index.html`, `css/app.css` | Page markup and styles |
| `js/gas-data.js` | Gas database, equation of state, mixture rules |
| `js/shock-calc.js` | Shock tube relations (no DOM access) |
| `js/riemann-solver.js`, `js/rmi-growth.js` | Exact Riemann solver, RMI growth models |
| `js/euler-solver.js`, `js/time-integrators.js` | 1D Euler solver |
| `js/ui-common.js`, `js/app-*.js`, `js/xt-*.js` | User interface |
| `test/` | Physics tests against exact solutions |

It is a static site with no build step: open `index.html`, or serve the folder with `python3 -m http.server`.

## Tests

```
npm test
```

Requires Node 20 or newer; there are no dependencies to install.
