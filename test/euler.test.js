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
