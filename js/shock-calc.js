/**
 * Zero-dimensional shock tube relations (incident shock, reflected shock,
 * driver pressure), generalized to the stiffened-gas equation of state.
 *
 * For a stiffened gas every ideal-gas jump relation holds with the pressure p
 * replaced by p + pinf, so the same expressions cover gases and liquids.
 *
 * Pure functions only: no DOM access.
 */

const SHOCK_ATM = 101325;       // Pa
const SHOCK_PSI_PER_ATM = 14.696;
const SHOCK_RU = 8314.51;       // J/(kmol K)
const HIGH_TEMPERATURE_WARNING_K = 1000;
const EXTREME_DRIVER_PRESSURE_PA = 2e8;

class ShockInputError extends Error {
    constructor(field, message) {
        super(message);
        this.name = 'ShockInputError';
        this.field = field;
    }
}

function shockRequirePositive(value, field, label) {
    if (typeof value !== 'number' || !isFinite(value) || value <= 0) {
        throw new ShockInputError(field, label + ' must be a positive number.');
    }
}

function shockValidateGas(gas, field, label) {
    if (!gas || !isFinite(gas.gamma) || gas.gamma <= 1 || !isFinite(gas.mw) || gas.mw <= 0) {
        throw new ShockInputError(field, 'Select a valid ' + label + '.');
    }
}

/** Particle velocity behind a shock of Mach number m moving into gas at rest. */
function shockParticleVelocity(m, gamma, a1) {
    return 2 * a1 / (gamma + 1) * (m - 1 / m);
}

/**
 * Largest shock Mach number a driver can produce (infinite driver pressure).
 * The expansion can accelerate driver gas to at most 2*a4/(gamma4 - 1).
 */
function maxMachForDriver(driven, driver, t1, t4) {
    const a1 = Math.sqrt(driven.gamma * SHOCK_RU / driven.mw * t1);
    const a4 = Math.sqrt(driver.gamma * SHOCK_RU / driver.mw * t4);
    const k = (driven.gamma + 1) * a4 / ((driver.gamma - 1) * a1);
    return 0.5 * (k + Math.sqrt(k * k + 4));
}

/**
 * Compute all shock tube states for a desired incident shock Mach number.
 *
 * @param {Object} input
 * @param {number} input.mach - incident shock Mach number (> 1)
 * @param {number} input.p1 - driven pressure, Pa
 * @param {number} input.t1 - driven temperature, K
 * @param {number} input.t4 - driver temperature, K
 * @param {{gamma:number, mw:number, pinf?:number}} input.driven
 * @param {{gamma:number, mw:number, pinf?:number}} input.driver
 * @param {number} [input.areaFactor=1] - empirical area-contraction factor
 * @returns {Object} states 1, 2, 4, 5 plus wave speeds and warnings.
 *   When the Mach number cannot be reached, driver pressure and density are
 *   null and a warning with code 'unreachable' is included; all other states
 *   are still returned.
 * @throws {ShockInputError} for non-physical input
 */
function computeShockTube(input) {
    const m = input.mach, p1 = input.p1, t1 = input.t1, t4 = input.t4;
    const areaFactor = input.areaFactor === undefined ? 1 : input.areaFactor;

    if (typeof m !== 'number' || !isFinite(m)) throw new ShockInputError('mach', 'Mach number must be a number.');
    if (m <= 1) throw new ShockInputError('mach', 'Mach number must be greater than 1 for a shock wave.');
    shockRequirePositive(p1, 'p1', 'Driven pressure');
    shockRequirePositive(t1, 't1', 'Driven temperature');
    shockRequirePositive(t4, 't4', 'Driver temperature');
    shockRequirePositive(areaFactor, 'areaFactor', 'Area factor');
    shockValidateGas(input.driven, 'driven', 'driven gas');
    shockValidateGas(input.driver, 'driver', 'driver gas');

    const g = input.driven.gamma, pinf = input.driven.pinf || 0;
    const g4 = input.driver.gamma, pinf4 = input.driver.pinf || 0;
    const gp = g + 1, gm = g - 1;
    const r = SHOCK_RU / input.driven.mw;
    const r4 = SHOCK_RU / input.driver.mw;

    // State 1 (stiffened pressures are written pt = p + pinf)
    const pt1 = p1 + pinf;
    const a1 = Math.sqrt(g * r * t1);
    const d1 = pt1 / (r * t1);
    const w = a1 * m;

    // State 2: incident shock
    const pr21 = 1 + 2 * g / gp * (m * m - 1);
    const pt2 = pr21 * pt1;
    const dr21 = (1 + gp / gm * pr21) / (gp / gm + pr21);
    const d2 = dr21 * d1;
    const t2 = t1 * pr21 / dr21;
    const a2 = Math.sqrt(g * r * t2);
    const u2 = shockParticleVelocity(m, g, a1);

    // State 5: shock reflected from the end wall
    const pr52 = ((3 * g - 1) * pr21 - gm) / (gm * pr21 + gp);
    const pt5 = pr52 * pt2;
    const dr52 = (1 + gp / gm * pr52) / (gp / gm + pr52);
    const d5 = dr52 * d2;
    const t5 = t2 * pr52 / dr52;
    const a5 = Math.sqrt(g * r * t5);
    const wr = d2 * u2 / (d5 - d2);
    const mr = (wr + u2) / a2;

    // State 4: isentropic expansion of the driver from rest to (u2, p2)
    const a4 = Math.sqrt(g4 * r4 * t4);
    const p2 = pt2 - pinf;
    const expansionArg = 1 - (g4 - 1) * u2 / (2 * a4);
    const warnings = [];
    let p4 = null, d4 = null, p4psig = null;

    if (expansionArg > 0) {
        p4 = (p2 + pinf4) * Math.pow(expansionArg, -2 * g4 / (g4 - 1)) - pinf4;
        if (!isFinite(p4)) p4 = null;
    }
    if (p4 === null) {
        const mMax = maxMachForDriver(input.driven, input.driver, t1, t4);
        warnings.push({
            code: 'unreachable',
            message: 'Mach ' + m + ' is unreachable with this driver (limit Mach ' + mMax.toFixed(2) +
                '). Use a lighter or hotter driver gas. Other states are still valid.'
        });
    } else {
        d4 = (p4 + pinf4) / (r4 * t4);
        p4psig = (p4 - SHOCK_ATM) / SHOCK_ATM * SHOCK_PSI_PER_ATM / areaFactor;
        if (p4 > EXTREME_DRIVER_PRESSURE_PA && pinf4 === 0) {
            warnings.push({
                code: 'extreme-driver-pressure',
                message: 'Driver pressure is impractically high (driver limit is Mach ' +
                    maxMachForDriver(input.driven, input.driver, t1, t4).toFixed(2) + '). Use a lighter or hotter driver gas.'
            });
        }
    }

    const hottest = Math.max(t2, t5);
    if (hottest > HIGH_TEMPERATURE_WARNING_K && pinf === 0) {
        warnings.push({
            code: 'high-temperature',
            message: 'Temperatures reach ' + hottest.toFixed(0) + ' K; constant γ overestimates them above about ' +
                HIGH_TEMPERATURE_WARNING_K + ' K.'
        });
    }

    return {
        mach: m,
        gamma: g, mw: input.driven.mw, pinf: pinf,
        state1: { p: p1, T: t1, rho: d1, a: a1, u: 0 },
        state2: { p: p2, T: t2, rho: d2, a: a2, u: u2 },
        state4: { p: p4, T: t4, rho: d4, a: a4, u: 0 },
        state5: { p: pt5 - pinf, T: t5, rho: d5, a: a5, u: 0 },
        shockSpeed: w,
        reflectedShockSpeed: wr,
        reflectedMach: mr,
        driverPressurePsig: p4psig,
        areaFactor: areaFactor,
        reachable: p4 !== null,
        warnings: warnings
    };
}

/**
 * Inverse problem: shock Mach number produced by a given driver pressure.
 * @param {Object} input - as computeShockTube, with p4 (Pa) instead of mach
 * @returns {number} Mach number
 * @throws {ShockInputError} if p4 does not exceed p1
 */
function machFromDriverPressure(input) {
    shockRequirePositive(input.p4, 'p4', 'Driver pressure');
    shockRequirePositive(input.p1, 'p1', 'Driven pressure');
    shockRequirePositive(input.t1, 't1', 'Driven temperature');
    shockRequirePositive(input.t4, 't4', 'Driver temperature');
    shockValidateGas(input.driven, 'driven', 'driven gas');
    shockValidateGas(input.driver, 'driver', 'driver gas');
    if (input.p4 <= input.p1) {
        throw new ShockInputError('p4', 'Driver pressure must exceed the driven pressure to produce a shock.');
    }

    const mMax = maxMachForDriver(input.driven, input.driver, input.t1, input.t4);
    const residual = function (m) {
        const res = computeShockTube({
            mach: m, p1: input.p1, t1: input.t1, t4: input.t4,
            driven: input.driven, driver: input.driver
        });
        return res.state4.p === null ? Infinity : res.state4.p - input.p4;
    };

    // p4(M) increases monotonically from p1 at M = 1 to infinity at mMax
    let lo = 1, hi = mMax;
    for (let i = 0; i < 200 && hi - lo > 1e-13 * hi; i++) {
        const mid = 0.5 * (lo + hi);
        if (mid === lo || mid === hi) break;
        if (residual(mid) > 0) hi = mid; else lo = mid;
    }
    return 0.5 * (lo + hi);
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        ShockInputError, computeShockTube, machFromDriverPressure,
        maxMachForDriver, shockParticleVelocity
    };
}
