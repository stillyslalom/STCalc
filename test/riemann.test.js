const test = require('node:test');
const assert = require('node:assert/strict');
const lib = require('./load.js');

const AIR = lib.GAS_DATABASE['1'];
const SF6 = lib.GAS_DATABASE['8'];
const WATER = lib.GAS_DATABASE['water'];

function close(actual, expected, rel, label) {
    assert.ok(Math.abs(actual - expected) <= rel * Math.abs(expected),
        `${label || 'value'}: expected ${expected}, got ${actual}`);
}

test('Sod shock tube (Toro test 1)', () => {
    // rho = 1 / 0.125, p = 1 / 0.1, gamma = 1.4 in non-dimensional units
    const gas = { gamma: 1.4, mw: lib.RU };          // R = 1
    const sol = lib.solveExactRiemann(
        { ...gas, p: 1, T: 1, u: 0 },
        { ...gas, p: 0.1, T: 0.8, u: 0 });
    close(sol.pStar, 0.30313, 1e-4, 'p*');
    close(sol.uStar, 0.92745, 1e-4, 'u*');
    close(sol.leftPost.rho, 0.42632, 1e-4, 'rho*L');
    close(sol.rightPost.rho, 0.26557, 1e-4, 'rho*R');
    assert.deepEqual({ ...sol.waveConfig }, { left: 'rarefaction', right: 'shock' });
});

test('symmetric states give a stationary contact', () => {
    const s = { ...AIR, p: 2e5, T: 350, u: 0 };
    const sol = lib.solveExactRiemann(s, s);
    close(sol.pStar, 2e5, 1e-9);
    assert.ok(Math.abs(sol.uStar) < 1e-6);
});

test('shocked air hitting SF6: reflected and transmitted shocks', () => {
    const shock = lib.computeShockTube({ mach: 1.5, p1: 101325, t1: 300, t4: 300, driven: AIR, driver: AIR });
    const sol = lib.solveShockedInterface(
        { ...AIR, p: shock.state2.p, T: shock.state2.T, u: shock.state2.u },
        { ...SF6, p: 101325, T: 300, u: 0 });
    assert.equal(sol.waveConfig.left, 'shock');
    assert.equal(sol.waveConfig.right, 'shock');
    assert.ok(sol.interfacePressure > shock.state2.p);
    assert.ok(sol.interfaceVelocity > 0 && sol.interfaceVelocity < shock.state2.u);
    assert.ok(sol.atwoodNumber > 0.6 && sol.atwoodNumber < 0.8);
});

test('air shock on water: nearly rigid reflection', () => {
    const shock = lib.computeShockTube({ mach: 2, p1: 101325, t1: 293.15, t4: 300, driven: AIR, driver: lib.GAS_DATABASE['4'] });
    const sol = lib.solveShockedInterface(
        { ...AIR, p: shock.state2.p, T: shock.state2.T, u: shock.state2.u },
        { ...WATER, p: 101325, T: 293.15, u: 0 });
    assert.equal(sol.converged, true);
    // Interface pressure approaches the rigid-wall reflected pressure
    close(sol.interfacePressure, shock.state5.p, 0.01, 'p interface');
    assert.ok(sol.interfaceVelocity > 0 && sol.interfaceVelocity < 1, `u = ${sol.interfaceVelocity}`);
    // Transmitted wave obeys the acoustic relation dp = rho * a * du in water
    const z = 998.2 * lib.eosSoundSpeed(WATER, 293.15);
    close(sol.interfacePressure - 101325, z * sol.interfaceVelocity, 0.01, 'acoustic impedance');
    close(sol.heavyGasPost.rho, 998.2, 0.001, 'water density');
});

test('high-pressure water expanding into air', () => {
    const sol = lib.solveExactRiemann(
        { ...WATER, p: 1e9, T: 293.15, u: 0 },
        { ...AIR, p: 1e5, T: 293.15, u: 0 });
    assert.equal(sol.converged, true);
    assert.deepEqual({ ...sol.waveConfig }, { left: 'rarefaction', right: 'shock' });
    assert.ok(sol.pStar > 1e5 && sol.pStar < 1e9);
    assert.ok(sol.uStar > 0);
});
