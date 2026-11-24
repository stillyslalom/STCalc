// Loads the browser scripts (plain globals, no module system) into one
// sandbox so the physics can be tested under Node without a DOM.
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const FILES = [
    '../data/nasa9-species.js',
    'nasa9-data.js',
    'thermo.js',
    'gas-data.js',
    'shock-calc.js',
    'riemann-solver.js',
    'time-integrators.js',
    'euler-solver.js',
    'xt-analysis.js'
];

function load() {
    const quiet = { log() {}, warn() {}, error() {} };
    const context = vm.createContext({ console: quiet, Math, Date });
    for (const file of FILES) {
        const source = fs.readFileSync(path.join(__dirname, '..', 'js', file), 'utf8');
        vm.runInContext(source, context, { filename: file });
    }
    return vm.runInContext(`({
        RU, GAS_DATABASE, stiffenedEffectiveMw, eosDensity, eosSoundSpeed,
        mixtureProperties, validateGasDefinition, gasComposition,
        NASA9_DATA, NASA9_ALL_SPECIES, NASA_ALIASES, searchNasaSpecies, nasa9Species, makeThermo,
        ShockInputError, computeShockTube, machFromDriverPressure, maxMachForDriver, tailoredMach,
        solveExactRiemann, sampleExactRiemann, solveShockedInterface,
        EulerSolver, steadyWindow, interfaceArrival, analyzeTestTime
    })`, context);
}

module.exports = load();
