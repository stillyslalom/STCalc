const test = require('node:test');
const assert = require('node:assert/strict');
const lib = require('./load.js');

const AIR = lib.GAS_DATABASE['1'];
const HELIUM = lib.GAS_DATABASE['4'];
const WATER = lib.GAS_DATABASE['water'];

function close(actual, expected, rel, label) {
    assert.ok(Math.abs(actual - expected) <= rel * Math.abs(expected),
        `${label || 'value'}: expected ${expected}, got ${actual}`);
}

const base = { p1: 101325, t1: 300, t4: 300, driven: AIR, driver: AIR, areaFactor: 1.35 };

test('air/air Mach 1.53 reproduces the original calculator', () => {
    const r = lib.computeShockTube({ ...base, mach: 1.53 });
    // Reference values are the page defaults of the original WiSTL calculator
    close(r.state2.p, 0.2599304145587011e6, 1e-9, 'p2');
    close(r.state2.T, 402.5381757242417, 1e-9, 'T2');
    close(r.state2.rho, 2.24958232829663, 1e-9, 'rho2');
    close(r.state2.u, 253.55410935058, 1e-9, 'u2');
    close(r.state4.p, 0.7858588290526163e6, 1e-9, 'p4');
    close(r.state4.rho, 9.125888293412393, 1e-9, 'rho4');
    close(r.state5.p, 0.5922303848754988e6, 1e-9, 'p5');
    close(r.state5.T, 517.9845985458894, 1e-9, 'T5');
    close(r.shockSpeed, 531.6185564241225, 1e-9, 'W');
    close(r.reflectedShockSpeed, 329.02882305300926, 1e-9, 'Wr');
    close(r.reflectedMach, 1.447458585280931, 1e-9, 'Mr');
    close(r.driverPressurePsig, 73.54339557717464, 1e-9, 'p4 psig');
    assert.equal(r.reachable, true);
    assert.equal(r.warnings.length, 0);
});

test('incident shock satisfies the Rankine-Hugoniot conditions', () => {
    for (const driven of [AIR, WATER]) {
        const r = lib.computeShockTube({ ...base, driven, driver: HELIUM, mach: 1.2 });
        const s1 = r.state1, s2 = r.state2, w = r.shockSpeed;
        const pinf = driven.pinf || 0;
        // Shock-fixed frame: gas enters at w, leaves at w - u2
        close(s1.rho * w, s2.rho * (w - s2.u), 1e-10, 'mass');
        close(s1.p + s1.rho * w * w, s2.p + s2.rho * (w - s2.u) ** 2, 1e-10, 'momentum');
        const h = (s, v) => driven.gamma * (s.p + pinf) / ((driven.gamma - 1) * s.rho) + 0.5 * v * v;
        close(h(s1, w), h(s2, w - s2.u), 1e-10, 'energy');
    }
});

test('reflected shock brings the gas to rest', () => {
    const r = lib.computeShockTube({ ...base, mach: 2.5, driver: HELIUM });
    const s2 = r.state2, s5 = r.state5, wr = r.reflectedShockSpeed;
    close(s2.rho * (s2.u + wr), s5.rho * wr, 1e-10, 'mass');
    close(s2.p + s2.rho * (s2.u + wr) ** 2, s5.p + s5.rho * wr * wr, 1e-10, 'momentum');
});

test('subsonic and non-numeric Mach numbers are rejected', () => {
    for (const mach of [0.5, 1, 0, -2, NaN]) {
        assert.throws(() => lib.computeShockTube({ ...base, mach }),
            err => err instanceof lib.ShockInputError && err.field === 'mach');
    }
    assert.throws(() => lib.computeShockTube({ ...base, mach: 2, p1: -1 }), e => e.field === 'p1');
    assert.throws(() => lib.computeShockTube({ ...base, mach: 2, t4: 0 }), e => e.field === 't4');
});

test('unreachable Mach number warns instead of failing', () => {
    const r = lib.computeShockTube({ ...base, mach: 8 });
    assert.equal(r.reachable, false);
    assert.equal(r.state4.p, null);
    assert.equal(r.driverPressurePsig, null);
    assert.ok(r.warnings.some(w => w.code === 'unreachable'));
    // The remaining states are still computed
    assert.ok(isFinite(r.state2.p) && isFinite(r.state5.T));
    const mMax = lib.maxMachForDriver(AIR, AIR, 300, 300);
    assert.ok(mMax > 6 && mMax < 6.5, `air/air limit ${mMax}`);
    assert.equal(lib.computeShockTube({ ...base, mach: mMax - 0.01 }).reachable, true);
});

test('absurd driver pressures are flagged', () => {
    const r = lib.computeShockTube({ ...base, mach: 5 });
    assert.equal(r.reachable, true);
    assert.ok(r.warnings.some(w => w.code === 'extreme-driver-pressure'));
});

test('high temperatures raise a constant-gamma warning', () => {
    const r = lib.computeShockTube({ ...base, mach: 4, driver: HELIUM });
    assert.ok(r.warnings.some(w => w.code === 'high-temperature'));
});

test('inverse solve recovers the Mach number', () => {
    for (const [driver, mach] of [[AIR, 1.53], [HELIUM, 3.2], [HELIUM, 1.01]]) {
        const forward = lib.computeShockTube({ ...base, driver, mach });
        const back = lib.machFromDriverPressure({ ...base, driver, p4: forward.state4.p });
        close(back, mach, 1e-8, 'Mach');
    }
    assert.throws(() => lib.machFromDriverPressure({ ...base, p4: 50000 }), e => e.field === 'p4');
});

test('water preset matches its reference state', () => {
    close(lib.eosDensity(WATER, 101325, 293.15), 998.2, 1e-12, 'density');
    const a = lib.eosSoundSpeed(WATER, 293.15);
    assert.ok(a > 1400 && a < 1500, `sound speed ${a}`);
});

test('mixture properties', () => {
    const gases = lib.GAS_DATABASE;
    const mix = lib.mixtureProperties([{ gasId: '4', fraction: 0.5 }, { gasId: '2', fraction: '0.5' }], gases);
    close(mix.mw, 0.5 * (4.002602 + 39.948), 1e-12, 'mw');
    assert.ok(mix.gamma > 1.667 && mix.gamma < 1.67);
    assert.throws(() => lib.mixtureProperties([{ gasId: '1', fraction: 0.4 }], gases), /sum to 1/);
    assert.throws(() => lib.mixtureProperties([], gases), /No gas components/);
    assert.throws(() => lib.mixtureProperties(
        [{ gasId: '1', fraction: 0.5 }, { gasId: 'water', fraction: 0.5 }], gases), /cannot be blended/);
    assert.equal(lib.mixtureProperties([{ gasId: 'water', fraction: 1 }], gases).pinf, 3.43e8);
});

/* ---- Area contraction at the diaphragm ---- */

test('area ratio 1 is the constant-area tube and larger ratios need less driver pressure', () => {
    const input = { ...base, driver: HELIUM, mach: 3 };
    const constant = lib.computeShockTube(input).state4.p;
    close(lib.computeShockTube({ ...input, areaRatio: 1.0001 }).state4.p, constant, 2e-4, 'continuity at 1');
    let previous = constant;
    for (const areaRatio of [1.2, 1.5, 2.25, 4, 10, 100]) {
        const p4 = lib.computeShockTube({ ...input, areaRatio }).state4.p;
        assert.ok(p4 < previous, `p4 falls with area ratio ${areaRatio}`);
        previous = p4;
    }
    // States outside the driver do not depend on the area ratio
    close(lib.computeShockTube({ ...input, areaRatio: 4 }).state5.p, lib.computeShockTube(input).state5.p, 1e-12, 'p5');
});

test('very large area ratio matches the closed-form steady-then-unsteady expansion', () => {
    const g1 = AIR.gamma, g4 = HELIUM.gamma;
    const a1 = Math.sqrt(g1 * lib.RU / AIR.mw * 300), a4 = Math.sqrt(g4 * lib.RU / HELIUM.mw * 300);
    for (const mach of [2, 3, 5]) {
        const u2 = 2 * a1 / (g1 + 1) * (mach - 1 / mach);
        const p2 = 101325 * (1 + 2 * g1 / (g1 + 1) * (mach * mach - 1));
        const aStar = a4 * Math.sqrt(2 / (g4 + 1));
        let a3;
        if (u2 > aStar) {
            // Steady expansion from rest to Mach 1, then unsteady from the sonic state to u2
            a3 = aStar * (g4 + 1) / 2 - (g4 - 1) / 2 * u2;
        } else {
            // Steady expansion from rest all the way to u2
            a3 = Math.sqrt(a4 * a4 - (g4 - 1) / 2 * u2 * u2);
        }
        const expected = p2 * Math.pow(a4 / a3, 2 * g4 / (g4 - 1));
        close(lib.computeShockTube({ ...base, driver: HELIUM, mach, areaRatio: 1e7 }).state4.p, expected, 1e-5, `Mach ${mach}`);
    }
});

test('area ratio with a subsonic nozzle exit conserves mass and energy', () => {
    // Weak shock: the nozzle does not choke. Check the steady-flow relations directly.
    const mach = 1.15, areaRatio = 3;
    const r = lib.computeShockTube({ ...base, driver: AIR, mach, areaRatio });
    const g = AIR.gamma, R = lib.RU / AIR.mw, a4 = Math.sqrt(g * R * 300);
    const T3 = 300 * Math.pow(r.state2.p / r.state4.p, (g - 1) / g);       // isentropic
    const a3 = Math.sqrt(g * R * T3), u3 = r.state2.u;
    assert.ok(u3 < a3, 'exit is subsonic');
    // Driver-side state from stagnation enthalpy and continuity
    const H0 = g * R * T3 / (g - 1) + 0.5 * u3 * u3;
    const rho = T => Math.pow(T, 1 / (g - 1));
    let lo = T3, hi = 300;
    for (let i = 0; i < 100; i++) {
        const T = 0.5 * (lo + hi);
        const ua = Math.sqrt(2 * (H0 - g * R * T / (g - 1)));
        if (areaRatio * rho(T) * ua > rho(T3) * u3) lo = T; else hi = T;
    }
    const Ta = 0.5 * (lo + hi), ua = Math.sqrt(2 * (H0 - g * R * Ta / (g - 1)));
    // That state must lie on the unsteady expansion from rest: ua = 2 (a4 - aa) / (gamma - 1)
    close(ua, 2 * (a4 - Math.sqrt(g * R * Ta)) / (g - 1), 1e-6, 'unsteady expansion relation');
});

test('area ratio: gas models agree for a monatomic driver, and the inverse solve works', () => {
    const helium = { gamma: 5 / 3, mw: 4.002602, thermo: [{ species: 'He', x: 1 }] };
    for (const areaRatio of [1.2, 2.25, 10]) {
        for (const mach of [1.3, 3, 6]) {
            const input = { ...base, driver: helium, mach, areaRatio };
            close(lib.computeShockTube({ ...input, variableGamma: true }).state4.p,
                lib.computeShockTube(input).state4.p, 1e-8, `area ratio ${areaRatio}, Mach ${mach}`);
        }
    }
    const forward = lib.computeShockTube({ ...base, driver: HELIUM, mach: 2.6, areaRatio: 2.25 });
    close(lib.machFromDriverPressure({ ...base, driver: HELIUM, p4: forward.state4.p, areaRatio: 2.25 }), 2.6, 1e-7, 'Mach');
    assert.ok(lib.maxMachForDriver(AIR, AIR, 300, 300, false, 4) > lib.maxMachForDriver(AIR, AIR, 300, 300, false, 1));
    assert.throws(() => lib.computeShockTube({ ...base, mach: 2, areaRatio: 0.5 }), e => e.field === 'areaRatio');
});
