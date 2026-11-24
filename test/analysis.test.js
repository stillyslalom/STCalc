const test = require('node:test');
const assert = require('node:assert/strict');
const lib = require('./load.js');

const AIR = lib.GAS_DATABASE['1'], HELIUM = lib.GAS_DATABASE['4'];

function close(actual, expected, rel, label) {
    assert.ok(Math.abs(actual - expected) <= rel * Math.abs(expected),
        `${label || 'value'}: expected ${expected}, got ${actual}`);
}

test('steadyWindow finds the arrival, plateau and end of a step', () => {
    const times = [], values = [];
    for (let j = 0; j <= 100; j++) {
        const t = j * 0.1;
        times.push(t);
        // Rise spread over three samples from t = 2, plateau of 5, second wave at t = 6
        values.push(t < 2 ? 1 : t < 2.3 ? 1 + (t - 2) / 0.3 * 4 : t < 6 ? 5 : 8);
    }
    const w = lib.steadyWindow(times, values);
    close(w.arrival, 2.1, 1e-9, 'arrival');
    close(w.plateau, 5, 1e-9, 'plateau');
    close(w.start + w.duration, 6, 1e-9, 'end');
    assert.equal(w.open, false);

    assert.equal(lib.steadyWindow(times, times.map(() => 3)), null);
    const open = lib.steadyWindow(times, times.map(t => t < 4 ? 1 : 2));
    assert.equal(open.open, true);
    close(open.duration, 6, 1e-9, 'open duration');
});

test('interfaceArrival', () => {
    const trajectory = [{ t: 0, x: 3 }, { t: 1, x: 4 }, { t: 2, x: 5.5 }, { t: 3, x: 5.2 }];
    assert.equal(lib.interfaceArrival(trajectory, 5), 2);
    assert.equal(lib.interfaceArrival(trajectory, 7), null);
});

test('tailoring: the reflected shock crosses a tailored interface without a reflected wave', () => {
    const base = { p1: 101325, t1: 300, t4: 300, driven: AIR, driver: HELIUM };
    const mach = lib.tailoredMach(base);
    assert.ok(mach > 3 && mach < 4, `tailored Mach ${mach}`);
    const check = m => {
        const r = lib.computeShockTube({ ...base, mach: m });
        // Stopped, reshocked driven gas (right) meets expanded driver gas still moving at u2 (left)
        const riemann = lib.solveExactRiemann(
            { ...HELIUM, p: r.state3.p, T: r.state3.T, u: r.state3.u },
            { ...AIR, p: r.state5.p, T: r.state5.T, u: 0 });
        return { r, riemann };
    };
    const tailored = check(mach);
    close(tailored.r.tailoringRatio, 1, 1e-9, 'ratio');
    close(tailored.riemann.pStar, tailored.r.state5.p, 1e-6, 'pressure is unchanged');
    assert.ok(Math.abs(tailored.riemann.uStar) < 1e-3, 'interface stops');

    const over = check(mach + 1), under = check(mach - 1);
    assert.ok(over.r.tailoringRatio > 1 && over.riemann.pStar > over.r.state5.p, 'over-tailored: pressure rises');
    assert.ok(under.r.tailoringRatio < 1 && under.riemann.pStar < under.r.state5.p, 'under-tailored: pressure falls');

    // State 3 is on the driver isentrope at the contact pressure and velocity
    const s3 = tailored.r.state3, s4 = tailored.r.state4;
    close(s3.p, tailored.r.state2.p, 1e-12, 'p3 = p2');
    close(s3.T, 300 * Math.pow(s3.p / s4.p, (HELIUM.gamma - 1) / HELIUM.gamma), 1e-9, 'T3');

    // The same gas on both sides is never tailored
    assert.equal(lib.tailoredMach({ ...base, driver: AIR }), null);
    // Both gas models agree for helium driving argon (both monatomic)
    const ar = lib.GAS_DATABASE['2'];
    const mono = id => ({ gamma: 5 / 3, mw: lib.GAS_DATABASE[id].mw, thermo: lib.gasComposition(lib.GAS_DATABASE[id]) });
    close(lib.tailoredMach({ ...base, driven: mono('2'), driver: mono('4'), variableGamma: true }),
        lib.tailoredMach({ ...base, driven: mono('2'), driver: mono('4') }), 1e-5, 'gas models');
    assert.ok(ar);
});

test('x-t end-wall test time: long at the tailored Mach number, short away from it', () => {
    const base = { p1: 101325, t1: 300, t4: 300, driven: AIR, driver: HELIUM };
    const slab = (gas, pressure, length) => ({ gasId: gas.name, gamma: gas.gamma, mw: gas.mw, pinf: 0, pressure, temperature: 300, length });
    const simulate = mach => {
        const shock = lib.computeShockTube({ ...base, mach });
        const solver = new lib.EulerSolver({
            nx: 600, cfl: 0.8, integrator: 'Hancock', order: 2, interfaceMethod: 'mixed',
            finalTime: 6 / shock.shockSpeed * 2.2, xtInterval: 2e-5
        });
        // Long driver so that the expansion reflected from the driver end wall arrives late
        solver.initialize([slab(HELIUM, shock.state4.p, 8), slab(AIR, 101325, 6)]);
        solver.run();
        const analysis = lib.analyzeTestTime(solver.getResults(), [13]);   // 5 m into the driven section
        return { shock, analysis, tEnd: solver.t };
    };

    const tailored = simulate(lib.tailoredMach(base));
    const wall = tailored.analysis.endWall;
    close(wall.arrival, 6 / tailored.shock.shockSpeed, 0.03, 'shock reaches the end wall');
    close(wall.plateau, tailored.shock.state5.p, 0.02, 'reflected pressure');

    const off = simulate(2.2);
    close(off.analysis.endWall.plateau, off.shock.state5.p, 0.02, 'reflected pressure, off-tailored');
    assert.equal(off.analysis.endWall.open, false, 'the interface interaction ends the test time');
    // Relative to the time the shock takes to cross the tube, tailoring gives a much longer test time
    const relative = s => s.analysis.endWall.duration / (6 / s.shock.shockSpeed);
    assert.ok(relative(tailored) > 2 * relative(off), `tailored ${relative(tailored)} vs ${relative(off)}`);

    // Probe: incident shock arrival and the state behind it
    const probe = tailored.analysis.probes[0];
    close(probe.window.arrival, 5 / tailored.shock.shockSpeed, 0.03, 'shock reaches the probe');
    close(probe.window.plateau, tailored.shock.state2.p, 0.02, 'post-shock pressure at the probe');
    assert.ok(probe.interfaceArrival === null || probe.interfaceArrival > probe.window.arrival);
});
