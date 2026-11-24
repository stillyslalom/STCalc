const test = require('node:test');
const assert = require('node:assert/strict');
const lib = require('./load.js');

const gas = id => ({ ...lib.GAS_DATABASE[id], thermo: lib.gasComposition(lib.GAS_DATABASE[id]) });
const AIR = gas('1'), ARGON = gas('2'), HELIUM = gas('4'), N2 = gas('6'), SF6 = gas('8');

function close(actual, expected, rel, label) {
    assert.ok(Math.abs(actual - expected) <= rel * Math.abs(expected),
        `${label || 'value'}: expected ${expected}, got ${actual}`);
}

test('NASA-9 properties at room temperature', () => {
    const n2 = lib.makeThermo(N2.thermo);
    close(n2.cp(298.15), 1039.7, 1e-3, 'cp of N2');        // J/(kg K)
    close(n2.gamma(300), 1.3995, 1e-3, 'gamma of N2');
    close(n2.a(300), 353.1, 2e-3, 'sound speed of N2');
    const air = lib.makeThermo(AIR.thermo);
    close(air.mw, 28.965, 1e-4, 'mw of air');
    close(air.gamma(300), 1.400, 1e-3, 'gamma of air');
    assert.ok(air.gamma(2000) < 1.30, 'gamma of air falls with temperature');
    close(lib.makeThermo(ARGON.thermo).gamma(2000), 5 / 3, 1e-8, 'gamma of argon');
});

test('enthalpy and entropy are consistent with cp', () => {
    const th = lib.makeThermo(AIR.thermo);
    for (const T of [150, 250, 600, 999, 1001, 3000, 7000]) {
        const dT = 1e-3 * T;
        close((th.h(T + dT) - th.h(T - dT)) / (2 * dT), th.cp(T), 1e-5, `dh/dT at ${T} K`);
        close((th.s0(T + dT) - th.s0(T - dT)) / (2 * dT), th.cp(T) / T, 1e-5, `ds/dT at ${T} K`);
        close(th.temperatureFromEnthalpy(th.h(T), 300), T, 1e-9, `inverse at ${T} K`);
    }
});

test('monatomic gases give the constant-gamma answer', () => {
    // Same gamma and molecular weight in both models
    const he = { gamma: 5 / 3, mw: 4.002602, thermo: HELIUM.thermo };
    const ar = { gamma: 5 / 3, mw: 39.948, thermo: ARGON.thermo };
    const input = { mach: 3.5, p1: 20000, t1: 295, t4: 295, driven: ar, driver: he, areaFactor: 1 };
    const a = lib.computeShockTube(input);
    const b = lib.computeShockTube({ ...input, variableGamma: true });
    assert.equal(b.variableGamma, true);
    for (const s of ['state2', 'state5']) {
        for (const k of ['p', 'T', 'rho', 'a']) close(b[s][k], a[s][k], 1e-8, `${s}.${k}`);
    }
    close(b.state2.u, a.state2.u, 1e-8, 'u2');
    close(b.reflectedShockSpeed, a.reflectedShockSpeed, 1e-8, 'Wr');
    close(b.state4.p, a.state4.p, 1e-6, 'p4');
    // Argon driving argon tops out near Mach 4.24, where argon is still calorically perfect
    close(lib.maxMachForDriver(ar, ar, 295, 295, true), lib.maxMachForDriver(ar, ar, 295, 295, false), 1e-6, 'Mach limit');
});

test('a user-defined constant-gamma gas is unchanged by the model switch', () => {
    const custom = { gamma: 1.31, mw: 16.04 };
    custom.thermo = lib.gasComposition(custom);
    const input = { mach: 2.4, p1: 101325, t1: 300, t4: 300, driven: custom, driver: custom, areaFactor: 1 };
    const a = lib.computeShockTube(input);
    const b = lib.computeShockTube({ ...input, variableGamma: true });
    close(b.state5.T, a.state5.T, 1e-8, 'T5');
    close(b.state4.p, a.state4.p, 1e-6, 'p4');
});

test('variable-gamma shocks satisfy the conservation laws', () => {
    for (const [driven, mach] of [[AIR, 1.5], [AIR, 6], [SF6, 3], [N2, 9]]) {
        const th = lib.makeThermo(driven.thermo);
        const r = lib.computeShockTube({
            mach, p1: 5000, t1: 300, t4: 300, driven, driver: HELIUM, variableGamma: true
        });
        const s1 = r.state1, s2 = r.state2, s5 = r.state5, w = r.shockSpeed, wr = r.reflectedShockSpeed;
        close(s1.rho * w, s2.rho * (w - s2.u), 1e-9, 'incident mass');
        close(s1.p + s1.rho * w * w, s2.p + s2.rho * (w - s2.u) ** 2, 1e-9, 'incident momentum');
        close(th.h(s2.T) - th.h(s1.T), 0.5 * (w * w - (w - s2.u) ** 2), 1e-8, 'incident energy');
        close(s2.rho * (s2.u + wr), s5.rho * wr, 1e-9, 'reflected mass');
        close(s2.p + s2.rho * (s2.u + wr) ** 2, s5.p + s5.rho * wr * wr, 1e-9, 'reflected momentum');
        close(th.h(s5.T) - th.h(s2.T), 0.5 * ((s2.u + wr) ** 2 - wr * wr), 1e-8, 'reflected energy');
        close(s2.p, s2.rho * th.R * s2.T, 1e-9, 'equation of state');
    }
});

test('hot air is cooler and denser than constant gamma predicts', () => {
    const input = { mach: 6, p1: 5000, t1: 300, t4: 300, driven: AIR, driver: HELIUM };
    const constant = lib.computeShockTube(input);
    const variable = lib.computeShockTube({ ...input, variableGamma: true });
    assert.ok(variable.state2.T < 0.9 * constant.state2.T, `T2 ${variable.state2.T} vs ${constant.state2.T}`);
    assert.ok(variable.state5.T < 0.85 * constant.state5.T, `T5 ${variable.state5.T} vs ${constant.state5.T}`);
    assert.ok(variable.state2.rho > constant.state2.rho);
    assert.ok(variable.state2.gamma < 1.33);
    // Weak shocks barely differ (the two models use slightly different gamma and mw)
    const weakC = lib.computeShockTube({ ...input, mach: 1.2 });
    const weakV = lib.computeShockTube({ ...input, mach: 1.2, variableGamma: true });
    close(weakV.state2.T, weakC.state2.T, 2e-3, 'weak-shock T2');
});

test('driver expansion with variable gamma is isentropic', () => {
    const r = lib.computeShockTube({
        mach: 2.5, p1: 101325, t1: 300, t4: 600, driven: AIR, driver: N2, variableGamma: true
    });
    // Compare against the closed form bracketed by gamma at the two end temperatures
    const lo = lib.computeShockTube({ mach: 2.5, p1: 101325, t1: 300, t4: 600, driven: AIR, driver: { gamma: 1.40, mw: 28.0134 } });
    const hi = lib.computeShockTube({ mach: 2.5, p1: 101325, t1: 300, t4: 600, driven: AIR, driver: { gamma: 1.37, mw: 28.0134 } });
    assert.ok(r.reachable && r.state4.p > 101325);
    const a = Math.min(lo.state4.p, hi.state4.p), b = Math.max(lo.state4.p, hi.state4.p);
    // Driven-side model differs between the runs, so allow a margin
    assert.ok(r.state4.p > 0.9 * a && r.state4.p < 1.1 * b, `p4 ${r.state4.p} outside [${a}, ${b}]`);
});

test('inverse solve and unreachable warning with variable gamma', () => {
    const base = { p1: 101325, t1: 300, t4: 300, driven: AIR, driver: HELIUM, variableGamma: true };
    const forward = lib.computeShockTube({ ...base, mach: 3 });
    close(lib.machFromDriverPressure({ ...base, p4: forward.state4.p }), 3, 1e-7, 'Mach');
    const far = lib.computeShockTube({ ...base, driver: AIR, mach: 9 });
    assert.equal(far.reachable, false);
    assert.equal(far.state4.p, null);
    assert.ok(far.warnings.some(w => w.code === 'unreachable'));
});
