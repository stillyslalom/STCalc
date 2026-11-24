const test = require('node:test');
const assert = require('node:assert/strict');
const lib = require('./load.js');

const AIR = lib.GAS_DATABASE['1'];
const HELIUM = lib.GAS_DATABASE['4'];
const WATER = lib.GAS_DATABASE['water'];

function slab(gas, pressure, temperature, length, velocity) {
    return { gasId: gas.name, gamma: gas.gamma, mw: gas.mw, pinf: gas.pinf || 0, pressure, temperature, length, velocity };
}

function run(slabs, options) {
    const solver = new lib.EulerSolver({
        nx: 400, cfl: 0.8, integrator: 'SSP', interfaceMethod: 'mixed', xtInterval: 1, ...options
    });
    solver.initialize(slabs);
    solver.run();
    return solver;
}

/** Mean of a field over cells whose center lies in [x0, x1]. */
function mean(solver, field, x0, x1) {
    let sum = 0, n = 0;
    for (let i = 0; i < solver.nx; i++) {
        if (solver.x[i] >= x0 && solver.x[i] <= x1) {
            sum += solver[field][i];
            n++;
        }
    }
    assert.ok(n > 0, 'no cells in range');
    return sum / n;
}

function close(actual, expected, rel, label) {
    assert.ok(Math.abs(actual - expected) <= rel * Math.abs(expected),
        `${label || 'value'}: expected ${expected}, got ${actual}`);
}

function totalMass(solver) {
    let m = 0;
    for (let i = 0; i < solver.nx; i++) m += solver.rho[i] * solver.dx;
    return m;
}

for (const method of ['mixed', 'sharp']) {
    for (const integrator of ['SSP', 'RK2']) {
        test(`helium/air shock tube matches the exact solution (${method}, ${integrator})`, () => {
            const left = { ...HELIUM, p: 4e5, T: 300, u: 0 };
            const right = { ...AIR, p: 101325, T: 300, u: 0 };
            const exact = lib.solveExactRiemann(left, right);

            const tEnd = 0.002;
            const solver = run([slab(HELIUM, 4e5, 300, 3), slab(AIR, 101325, 300, 3)],
                { finalTime: tEnd, interfaceMethod: method, integrator, cfl: integrator === 'SSP' ? 0.8 : 0.4 });

            // Sample the plateaus on either side of the contact
            const xContact = 3 + exact.uStar * tEnd;
            const shockSpeed = right.p < exact.pStar
                ? (exact.pStar - right.p) / (lib.eosDensity(AIR, right.p, right.T) * exact.uStar) : 0;
            const xShock = 3 + shockSpeed * tEnd;
            const a = xContact + 0.25 * (xShock - xContact), b = xContact + 0.75 * (xShock - xContact);
            // The sharp method is not conservative at the interface and is less accurate
            const tol = method === 'mixed' ? 0.01 : 0.03;
            close(mean(solver, 'p', a, b), exact.pStar, tol, 'p*');
            close(mean(solver, 'u', a, b), exact.uStar, tol, 'u*');
            close(mean(solver, 'rho', a, b), exact.rightPost.rho, 0.02, 'shocked air density');
            close(mean(solver, 'T', a, b), exact.rightPost.T, 0.02, 'shocked air temperature');
            close(solver.tracers[0].x, xContact, tol, 'contact position');
            close(totalMass(solver), totalMass(run([slab(HELIUM, 4e5, 300, 3), slab(AIR, 101325, 300, 3)],
                { finalTime: 1e-9, interfaceMethod: method })), 1e-10, 'mass conservation');
        });
    }
}

test('uniform state stays uniform', () => {
    const solver = run([slab(AIR, 101325, 300, 1), slab(AIR, 101325, 300, 1)], { finalTime: 0.001, nx: 100 });
    for (let i = 0; i < solver.nx; i++) {
        close(solver.p[i], 101325, 1e-12, 'p');
        assert.ok(Math.abs(solver.u[i]) < 1e-9);
    }
});

test('moving air/water interface keeps pressure and velocity uniform', () => {
    // Interface placed mid-cell so that the initial condition has a mixed cell
    const u0 = 50, p0 = 101325, tEnd = 0.0005;
    const solver = run([slab(AIR, p0, 293.15, 1.0025, u0), slab(WATER, p0, 293.15, 4.9975, u0)],
        { finalTime: tEnd, nx: 600 });
    // Short run and a window clear of the waves that the walls send into the moving fluid
    for (let i = 0; i < solver.nx; i++) {
        if (solver.x[i] < 1.0 || solver.x[i] > 1.3) continue;
        assert.ok(Math.abs(solver.p[i] - p0) < 1e-9 * p0, `p[${i}] = ${solver.p[i]}`);
        assert.ok(Math.abs(solver.u[i] - u0) < 1e-9 * u0, `u[${i}] = ${solver.u[i]}`);
    }
    close(solver.tracers[0].x, 1.0025 + u0 * tEnd, 1e-3, 'interface position');
    close(mean(solver, 'rho', 1.2, 1.3), 998.2, 1e-6, 'water density');
    close(mean(solver, 'T', 1.0, 1.05), 293.15, 1e-4, 'air temperature');
});

test('air shock reflecting from water matches the exact interface solution', () => {
    const shock = lib.computeShockTube({ mach: 2, p1: 101325, t1: 293.15, t4: 293.15, driven: AIR, driver: HELIUM });
    const exact = lib.solveShockedInterface(
        { ...AIR, p: shock.state2.p, T: shock.state2.T, u: shock.state2.u },
        { ...WATER, p: 101325, T: 293.15, u: 0 });

    // Shocked air (moving) | quiescent air | water. The shock forms at x = 1,
    // reaches the water at x = 1.5 and reflects.
    const tHit = 0.5 / shock.shockSpeed;
    const tEnd = tHit + 0.0004;
    const solver = run([
        slab(AIR, shock.state2.p, shock.state2.T, 1, shock.state2.u),
        slab(AIR, 101325, 293.15, 0.5),
        slab(WATER, 101325, 293.15, 1.5)
    ], { finalTime: tEnd, nx: 1200 });

    // Reflected shock has moved ~0.1 m back into the air; transmitted wave ~0.6 m into the water
    close(mean(solver, 'p', 1.44, 1.48), exact.interfacePressure, 0.03, 'air side pressure');
    close(mean(solver, 'p', 1.6, 1.9), exact.interfacePressure, 0.03, 'water side pressure');
    assert.ok(Math.abs(mean(solver, 'u', 1.6, 1.9) - exact.interfaceVelocity) < 0.1, 'water velocity');
    for (let i = 0; i < solver.nx; i++) {
        assert.ok(isFinite(solver.p[i]) && solver.rho[i] > 0, `cell ${i} is physical`);
    }
});

test('sharp interface method refuses gas next to liquid', () => {
    const solver = new lib.EulerSolver({ nx: 50, interfaceMethod: 'sharp' });
    assert.throws(() => solver.initialize([slab(AIR, 1e5, 300, 1), slab(WATER, 1e5, 300, 1)]), /Mixed-Cell/);
});

/* ---- Temperature-dependent (tabulated) equations of state ---- */

const comp = id => lib.gasComposition(lib.GAS_DATABASE[id]);
const tslab = (gas, thermo, pressure, temperature, length, velocity) =>
    ({ ...slab(gas, pressure, temperature, length, velocity), thermo });

test('tabulated constant-cp gases reproduce the constant-gamma solver', () => {
    const constant = id => [{ gamma: lib.GAS_DATABASE[id].gamma, mw: lib.GAS_DATABASE[id].mw, x: 1 }];
    const options = { finalTime: 0.002 };
    const a = run([slab(HELIUM, 4e5, 300, 3), slab(AIR, 101325, 300, 3)], options);
    const b = run([tslab(HELIUM, constant('4'), 4e5, 300, 3), tslab(AIR, constant('1'), 101325, 300, 3)], options);
    assert.equal(b.tabulated, true);
    for (let i = 0; i < a.nx; i++) {
        close(b.p[i], a.p[i], 1e-6, `p[${i}]`);
        close(b.rho[i], a.rho[i], 1e-6, `rho[${i}]`);
        close(b.T[i], a.T[i], 1e-6, `T[${i}]`);
        assert.ok(Math.abs(b.u[i] - a.u[i]) < 1e-4, `u[${i}]`);
    }
});

test('variable-gamma shock tube matches the 0D variable-gamma calculation', () => {
    const air = { ...AIR, thermo: comp('1') }, he = { ...HELIUM, thermo: comp('4') };
    const shock = lib.computeShockTube({ mach: 4, p1: 10000, t1: 300, t4: 300, driven: air, driver: he, variableGamma: true });
    const tEnd = 0.0012;
    const solver = run([tslab(HELIUM, he.thermo, shock.state4.p, 300, 3), tslab(AIR, air.thermo, 10000, 300, 3)],
        { finalTime: tEnd, nx: 800 });

    // Sample just behind the shock: the first-order scheme smears helium from
    // the contact across much of the short shocked region
    const xShock = 3 + shock.shockSpeed * tEnd;
    const a = xShock - 0.12, b = xShock - 0.06;
    close(mean(solver, 'p', a, b), shock.state2.p, 0.01, 'p2');
    close(mean(solver, 'u', a, b), shock.state2.u, 0.01, 'u2');
    close(mean(solver, 'T', a, b), shock.state2.T, 0.02, 'T2');
    close(mean(solver, 'rho', a, b), shock.state2.rho, 0.02, 'rho2');
    close(mean(solver, 'gamma', a, b), shock.state2.gamma, 0.01, 'gamma2');

    // Constant gamma gives a clearly hotter post-shock state for the same shock
    const constant = lib.computeShockTube({ mach: 4, p1: 10000, t1: 300, t4: 300, driven: AIR, driver: HELIUM });
    assert.ok(mean(solver, 'T', a, b) < 0.97 * constant.state2.T);
});

test('moving interface between variable-gamma air and water stays uniform', () => {
    const u0 = 50, p0 = 101325, tEnd = 0.0005;
    const solver = run([tslab(AIR, comp('1'), p0, 293.15, 1.0025, u0), slab(WATER, p0, 293.15, 4.9975, u0)],
        { finalTime: tEnd, nx: 600 });
    for (let i = 0; i < solver.nx; i++) {
        if (solver.x[i] < 1.0 || solver.x[i] > 1.3) continue;
        assert.ok(Math.abs(solver.p[i] - p0) < 1e-7 * p0, `p[${i}] = ${solver.p[i]}`);
        assert.ok(Math.abs(solver.u[i] - u0) < 1e-7 * u0, `u[${i}] = ${solver.u[i]}`);
    }
    close(mean(solver, 'T', 1.0, 1.02), 293.15, 1e-4, 'air temperature');
});

test('partial densities are conserved', () => {
    const solver = run([tslab(HELIUM, comp('4'), 8e5, 300, 2), tslab(AIR, comp('1'), 101325, 300, 2), tslab(lib.GAS_DATABASE['8'], comp('8'), 101325, 300, 1)],
        { finalTime: 0.004, nx: 300 });
    const M = solver.numMaterials;
    const expected = [8e5 / (lib.RU / 4.002602 * 300) * 2, null, null];
    let total = 0;
    const masses = new Array(M).fill(0);
    for (let i = 0; i < solver.nx; i++) {
        for (let m = 0; m < M; m++) masses[m] += solver.U[i * solver.nvar + 3 + M + m] * solver.dx;
        total += solver.rho[i] * solver.dx;
    }
    close(masses[0], expected[0], 1e-9, 'helium mass');
    close(masses[0] + masses[1] + masses[2], total, 1e-9, 'sum of partial masses');
    for (let i = 0; i < solver.nx; i++) assert.ok(isFinite(solver.p[i]) && solver.p[i] > 0 && solver.T[i] > 0);
});

test('sharp interface method refuses variable-gamma slabs', () => {
    const solver = new lib.EulerSolver({ nx: 50, interfaceMethod: 'sharp' });
    assert.throws(() => solver.initialize([tslab(AIR, comp('1'), 1e5, 300, 1), slab(AIR, 1e5, 300, 1)]), /Mixed-Cell/);
});

/* ---- Second-order (MUSCL) reconstruction ---- */

/** Number of cells across which density goes from 10% to 90% of a jump. */
function contactWidth(solver, rhoA, rhoB, x0, x1) {
    let cells = 0;
    for (let i = 0; i < solver.nx; i++) {
        const f = (solver.rho[i] - rhoA) / (rhoB - rhoA);
        if (solver.x[i] > x0 && solver.x[i] < x1 && f > 0.1 && f < 0.9) cells++;
    }
    return cells;
}

test('second order sharpens the contact and reduces the error', () => {
    const exact = lib.solveExactRiemann({ ...HELIUM, p: 4e5, T: 300, u: 0 }, { ...AIR, p: 101325, T: 300, u: 0 });
    const tEnd = 0.002, xContact = 3 + exact.uStar * tEnd;
    const slabs = [slab(HELIUM, 4e5, 300, 3), slab(AIR, 101325, 300, 3)];
    const first = run(slabs, { finalTime: tEnd, order: 1 });
    const second = run(slabs, { finalTime: tEnd, order: 2 });

    const width = s => contactWidth(s, exact.leftPost.rho, exact.rightPost.rho, xContact - 0.5, xContact + 0.3);
    assert.ok(width(second) <= 0.6 * width(first), `contact ${width(second)} vs ${width(first)} cells`);

    const i = Math.floor((xContact + 0.15) / second.dx);
    const error = s => Math.abs(s.u[i] / exact.uStar - 1);
    assert.ok(error(second) < 0.25 * error(first), `u* error ${error(second)} vs ${error(first)}`);
    close(second.p[i], exact.pStar, 1e-3, 'p*');
    close(second.tracers[0].x, xContact, 1e-3, 'contact position');
    close(totalMass(second), totalMass(first), 1e-10, 'mass');

    // No overshoot at the shock or the contact
    for (let k = 0; k < second.nx; k++) {
        assert.ok(second.p[k] <= 4e5 * (1 + 1e-9) && second.p[k] >= 101325 * (1 - 1e-9), `p[${k}] = ${second.p[k]}`);
    }
});

for (const variable of [false, true]) {
    test(`second order keeps a moving air/water interface uniform (${variable ? 'variable' : 'constant'} gamma)`, () => {
        const u0 = 50, p0 = 101325;
        const solver = run([tslab(AIR, variable ? comp('1') : null, p0, 293.15, 1.0025, u0), slab(WATER, p0, 293.15, 4.9975, u0)],
            { finalTime: 0.0005, nx: 600, order: 2 });
        for (let i = 0; i < solver.nx; i++) {
            if (solver.x[i] < 1.0 || solver.x[i] > 1.3) continue;
            assert.ok(Math.abs(solver.p[i] - p0) < 1e-7 * p0, `p[${i}] = ${solver.p[i]}`);
            assert.ok(Math.abs(solver.u[i] - u0) < 1e-7 * u0, `u[${i}] = ${solver.u[i]}`);
        }
    });
}

test('second order: air shock reflecting from water', () => {
    const shock = lib.computeShockTube({ mach: 2, p1: 101325, t1: 293.15, t4: 293.15, driven: AIR, driver: HELIUM });
    const exact = lib.solveShockedInterface(
        { ...AIR, p: shock.state2.p, T: shock.state2.T, u: shock.state2.u },
        { ...WATER, p: 101325, T: 293.15, u: 0 });
    const solver = run([
        slab(AIR, shock.state2.p, shock.state2.T, 1, shock.state2.u),
        slab(AIR, 101325, 293.15, 0.5),
        slab(WATER, 101325, 293.15, 1.5)
    ], { finalTime: 0.5 / shock.shockSpeed + 0.0004, nx: 600, order: 2 });
    close(mean(solver, 'p', 1.44, 1.48), exact.interfacePressure, 0.01, 'air side pressure');
    close(mean(solver, 'p', 1.6, 1.9), exact.interfacePressure, 0.01, 'water side pressure');
    for (let i = 0; i < solver.nx; i++) assert.ok(isFinite(solver.p[i]) && solver.rho[i] > 0);
});

test('second order with variable gamma matches the 0D calculation across the shocked region', () => {
    const air = { ...AIR, thermo: comp('1') }, he = { ...HELIUM, thermo: comp('4') };
    const shock = lib.computeShockTube({ mach: 4, p1: 10000, t1: 300, t4: 300, driven: air, driver: he, variableGamma: true });
    const tEnd = 0.0012;
    const solver = run([tslab(HELIUM, he.thermo, shock.state4.p, 300, 3), tslab(AIR, air.thermo, 10000, 300, 3)],
        { finalTime: tEnd, nx: 800, order: 2 });
    // The middle of the shocked region, which first order contaminates with helium
    const xShock = 3 + shock.shockSpeed * tEnd, xContact = 3 + shock.state2.u * tEnd;
    const a = xContact + 0.4 * (xShock - xContact), b = xContact + 0.8 * (xShock - xContact);
    close(mean(solver, 'T', a, b), shock.state2.T, 0.005, 'T2');
    close(mean(solver, 'p', a, b), shock.state2.p, 0.005, 'p2');
    close(mean(solver, 'rho', a, b), shock.state2.rho, 0.005, 'rho2');
});

test('second order works with the sharp interface method', () => {
    const solver = run([slab(HELIUM, 4e5, 300, 3), slab(AIR, 101325, 300, 3)],
        { finalTime: 0.002, order: 2, interfaceMethod: 'sharp' });
    const exact = lib.solveExactRiemann({ ...HELIUM, p: 4e5, T: 300, u: 0 }, { ...AIR, p: 101325, T: 300, u: 0 });
    close(mean(solver, 'p', 3.6, 3.9), exact.pStar, 0.03, 'p*');
});

/* ---- MUSCL-Hancock time stepping (one flux evaluation per step) ---- */

const hancock = { integrator: 'Hancock', order: 2, cfl: 0.8 };

test('Hancock: helium/air shock tube matches the exact solution', () => {
    const exact = lib.solveExactRiemann({ ...HELIUM, p: 4e5, T: 300, u: 0 }, { ...AIR, p: 101325, T: 300, u: 0 });
    const tEnd = 0.002, xContact = 3 + exact.uStar * tEnd;
    const slabs = [slab(HELIUM, 4e5, 300, 3), slab(AIR, 101325, 300, 3)];
    const solver = run(slabs, { finalTime: tEnd, ...hancock });
    const reference = run(slabs, { finalTime: tEnd, order: 2 });       // SSP Runge-Kutta

    const i = Math.floor((xContact + 0.15) / solver.dx);
    close(solver.p[i], exact.pStar, 1e-3, 'p*');
    close(solver.u[i], exact.uStar, 1e-3, 'u*');
    close(solver.rho[i], exact.rightPost.rho, 5e-3, 'shocked air density');
    close(solver.tracers[0].x, xContact, 1e-3, 'contact position');
    close(totalMass(solver), totalMass(reference), 1e-10, 'mass');
    const width = s => contactWidth(s, exact.leftPost.rho, exact.rightPost.rho, xContact - 0.5, xContact + 0.3);
    assert.ok(width(solver) <= width(reference) + 2, `contact ${width(solver)} vs ${width(reference)} cells`);
    for (let k = 0; k < solver.nx; k++) {
        assert.ok(solver.p[k] <= 4e5 * (1 + 1e-6) && solver.p[k] >= 101325 * (1 - 1e-6), `p[${k}] = ${solver.p[k]}`);
    }
    // One flux evaluation per step, at a similar number of steps
    assert.ok(solver.timeSteps < 1.1 * reference.timeSteps);
});

for (const variable of [false, true]) {
    test(`Hancock keeps a moving air/water interface uniform (${variable ? 'variable' : 'constant'} gamma)`, () => {
        const u0 = 50, p0 = 101325;
        const solver = run([tslab(AIR, variable ? comp('1') : null, p0, 293.15, 1.0025, u0), slab(WATER, p0, 293.15, 4.9975, u0)],
            { finalTime: 0.0005, nx: 600, ...hancock });
        for (let i = 0; i < solver.nx; i++) {
            if (solver.x[i] < 1.0 || solver.x[i] > 1.3) continue;
            assert.ok(Math.abs(solver.p[i] - p0) < 1e-7 * p0, `p[${i}] = ${solver.p[i]}`);
            assert.ok(Math.abs(solver.u[i] - u0) < 1e-7 * u0, `u[${i}] = ${solver.u[i]}`);
        }
        close(solver.tracers[0].x, 1.0025 + u0 * 0.0005, 1e-3, 'interface position');
    });
}

test('Hancock: air shock reflecting from water', () => {
    const shock = lib.computeShockTube({ mach: 2, p1: 101325, t1: 293.15, t4: 293.15, driven: AIR, driver: HELIUM });
    const exact = lib.solveShockedInterface(
        { ...AIR, p: shock.state2.p, T: shock.state2.T, u: shock.state2.u },
        { ...WATER, p: 101325, T: 293.15, u: 0 });
    const solver = run([
        slab(AIR, shock.state2.p, shock.state2.T, 1, shock.state2.u),
        slab(AIR, 101325, 293.15, 0.5),
        slab(WATER, 101325, 293.15, 1.5)
    ], { finalTime: 0.5 / shock.shockSpeed + 0.0004, nx: 600, ...hancock });
    close(mean(solver, 'p', 1.44, 1.48), exact.interfacePressure, 0.01, 'air side pressure');
    close(mean(solver, 'p', 1.6, 1.9), exact.interfacePressure, 0.01, 'water side pressure');
    for (let i = 0; i < solver.nx; i++) assert.ok(isFinite(solver.p[i]) && solver.rho[i] > 0);
});

test('Hancock with variable gamma matches the 0D calculation, including the reflected state', () => {
    const air = { ...AIR, thermo: comp('1') }, he = { ...HELIUM, thermo: comp('4') };
    const shock = lib.computeShockTube({ mach: 4, p1: 10000, t1: 300, t4: 300, driven: air, driver: he, variableGamma: true });
    // Shock reaches the end wall (3 m of driven gas) and reflects
    const tWall = 3 / shock.shockSpeed;
    const incident = run([tslab(HELIUM, he.thermo, shock.state4.p, 300, 3), tslab(AIR, air.thermo, 10000, 300, 3)],
        { finalTime: 0.6 * tWall, nx: 800, ...hancock });
    const xShock = 3 + shock.shockSpeed * 0.6 * tWall, xContact = 3 + shock.state2.u * 0.6 * tWall;
    const a = xContact + 0.4 * (xShock - xContact), b = xContact + 0.8 * (xShock - xContact);
    close(mean(incident, 'T', a, b), shock.state2.T, 0.005, 'T2');
    close(mean(incident, 'p', a, b), shock.state2.p, 0.005, 'p2');
    close(mean(incident, 'u', a, b), shock.state2.u, 0.005, 'u2');

    const reflected = run([tslab(HELIUM, he.thermo, shock.state4.p, 300, 3), tslab(AIR, air.thermo, 10000, 300, 3)],
        { finalTime: 1.04 * tWall, nx: 800, ...hancock });
    close(mean(reflected, 'p', 5.985, 5.995), shock.state5.p, 0.02, 'p5');
    close(mean(reflected, 'T', 5.985, 5.995), shock.state5.T, 0.03, 'T5');
});

test('Hancock with a first-order solver is Godunov\'s method', () => {
    const slabs = [slab(HELIUM, 4e5, 300, 3), slab(AIR, 101325, 300, 3)];
    const solver = run(slabs, { finalTime: 0.002, integrator: 'Hancock', order: 1, interfaceMethod: 'sharp' });
    const exact = lib.solveExactRiemann({ ...HELIUM, p: 4e5, T: 300, u: 0 }, { ...AIR, p: 101325, T: 300, u: 0 });
    close(mean(solver, 'p', 3.6, 3.9), exact.pStar, 0.03, 'p*');
    const mixed = run(slabs, { finalTime: 0.002, integrator: 'Hancock', order: 1 });
    close(mean(mixed, 'p', 3.6, 3.9), exact.pStar, 0.01, 'p* (mixed)');
});
