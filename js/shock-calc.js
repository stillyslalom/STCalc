/**
 * Zero-dimensional shock tube relations (incident shock, reflected shock,
 * driver pressure) for constant-gamma gases, stiffened-gas liquids and
 * thermally perfect gases with temperature-dependent specific heats.
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
const DISSOCIATION_WARNING_K = 2500;

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

/* ------------------------------------------------------------------ */
/* Gas models                                                          */
/*                                                                     */
/* Each side of the diaphragm is described by a "model" with the same  */
/* three operations, so that constant-gamma and temperature-dependent  */
/* gases can be mixed freely.                                          */
/* ------------------------------------------------------------------ */

/** Calorically perfect (or stiffened) gas: closed-form relations. */
function constantGammaModel(gas) {
    const g = gas.gamma, pinf = gas.pinf || 0;
    const gp = g + 1, gm = g - 1;
    const R = SHOCK_RU / gas.mw;
    const state = function (pt, T, u) {       // pt = p + pinf
        return { p: pt - pinf, T: T, rho: pt / (R * T), a: Math.sqrt(g * R * T), u: u, gamma: g };
    };
    const jump = function (ratio) {           // density ratio for a shifted-pressure ratio
        return (1 + gp / gm * ratio) / (gp / gm + ratio);
    };
    return {
        mw: gas.mw, pinf: pinf, thermo: null,
        state: function (p, T) { return state(p + pinf, T, 0); },
        /** State behind a shock of Mach number m moving into s1 at rest. */
        incidentShock: function (s1, m) {
            const ratio = 1 + 2 * g / gp * (m * m - 1);
            const dr = jump(ratio);
            return state(ratio * (s1.p + pinf), s1.T * ratio / dr, 2 * s1.a / gp * (m - 1 / m));
        },
        /** State behind the shock that brings s2 to rest at a wall, and its lab-frame speed. */
        reflectedShock: function (s1, s2) {
            const r21 = (s2.p + pinf) / (s1.p + pinf);
            const r52 = ((3 * g - 1) * r21 - gm) / (gm * r21 + gp);
            const dr = jump(r52);
            const s5 = state(r52 * (s2.p + pinf), s2.T * r52 / dr, 0);
            return { state: s5, speed: s2.rho * s2.u / (s5.rho - s2.rho) };
        },
        /** Fastest speed an unsteady expansion from rest at T4 can reach. */
        maxExpansionVelocity: function (T4) {
            return 2 * Math.sqrt(g * R * T4) / gm;
        },
        /** Initial pressure that expands isentropically from rest at T4 to (u, p). */
        expansionSourcePressure: function (T4, u, p) {
            const arg = 1 - gm * u / (2 * Math.sqrt(g * R * T4));
            if (!(arg > 0)) return null;
            const p4 = (p + pinf) * Math.pow(arg, -2 * g / gm) - pinf;
            return isFinite(p4) ? p4 : null;
        }
    };
}

/** Bisection for a root of f bracketed by [lo, hi]. */
function shockBisect(f, lo, hi) {
    let flo = f(lo);
    for (let i = 0; i < 200; i++) {
        const mid = 0.5 * (lo + hi);
        if (mid === lo || mid === hi) break;
        const fmid = f(mid);
        if ((fmid > 0) === (flo > 0)) { lo = mid; flo = fmid; } else { hi = mid; }
    }
    return 0.5 * (lo + hi);
}

/**
 * Thermally perfect ideal gas: cp depends on temperature (frozen composition).
 * Jump conditions are solved numerically for the density ratio.
 * @param {Object} th - property functions from makeThermo()
 */
function thermallyPerfectModel(th) {
    const R = th.R;
    const state = function (p, T, u) {
        return { p: p, T: T, rho: p / (R * T), a: th.a(T), u: u, gamma: th.gamma(T) };
    };

    /**
     * Normal shock with upstream state s and upstream speed v in the shock
     * frame. Mass, momentum and energy give p and h downstream for a trial
     * density ratio eps = rho_up / rho_down; the equation of state closes it.
     * eps = 1 (no shock) always satisfies the equations, so that root is
     * divided out.
     */
    const shock = function (s, v) {
        const h1 = th.h(s.T);
        const downstream = function (eps) {
            const p = s.p + s.rho * v * v * (1 - eps);
            const T = th.temperatureFromEnthalpy(h1 + 0.5 * v * v * (1 - eps * eps), s.T);
            return { p: p, T: T };
        };
        const residual = function (eps) {
            const d = downstream(eps);
            return (s.rho * R * d.T / d.p - eps) / (1 - eps);
        };
        const eps = shockBisect(residual, 1e-6, 1 - 1e-7);
        const d = downstream(eps);
        return { p: d.p, T: d.T, eps: eps };
    };

    // The expansion is integrated in sqrt(T), where the integrand of
    // u = integral of cp/a dT is smooth all the way down to T = 0.
    const expansionVelocity = function (T4, T3) {
        const lo = Math.sqrt(Math.max(T3, 0)), hi = Math.sqrt(T4);
        const n = 200, dx = (hi - lo) / n;
        const f = function (x) {
            const T = Math.max(x * x, 1e-12);
            return 2 * th.cp(T) / Math.sqrt(th.gamma(T) * R);
        };
        let sum = f(lo) + f(hi);
        for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * f(lo + i * dx);
        return sum * dx / 3;   // Simpson's rule
    };

    return {
        mw: th.mw, pinf: 0, thermo: th,
        state: function (p, T) { return state(p, T, 0); },
        incidentShock: function (s1, m) {
            const w = m * s1.a;
            const d = shock(s1, w);
            return state(d.p, d.T, w * (1 - d.eps));
        },
        reflectedShock: function (s1, s2) {
            // Shock frame: gas arrives at u2 + wr and leaves at wr, so with
            // eps = rho2/rho5 the upstream speed is u2 / (1 - eps).
            const h2 = th.h(s2.T);
            const downstream = function (eps) {
                const p = s2.p + s2.rho * s2.u * s2.u / (1 - eps);
                const T = th.temperatureFromEnthalpy(h2 + 0.5 * s2.u * s2.u * (1 + eps) / (1 - eps), s2.T);
                return { p: p, T: T };
            };
            const eps = shockBisect(function (e) {
                const d = downstream(e);
                return s2.rho * R * d.T / d.p - e;
            }, 1e-6, 1 - 1e-9);
            const d = downstream(eps);
            return { state: state(d.p, d.T, 0), speed: s2.u * eps / (1 - eps) };
        },
        maxExpansionVelocity: function (T4) {
            return expansionVelocity(T4, 0);
        },
        expansionSourcePressure: function (T4, u, p) {
            if (!(u < expansionVelocity(T4, 0))) return null;
            const T3 = shockBisect(function (T) { return expansionVelocity(T4, T) - u; }, 0, T4);
            // Isentropic: s0(T) - R ln(p) is constant
            const p4 = p * Math.exp((th.s0(T4) - th.s0(T3)) / R);
            return isFinite(p4) ? p4 : null;
        }
    };
}

/**
 * Model for one gas. The temperature-dependent model is used when requested
 * and the gas has a composition for it; liquids always use the stiffened gas.
 */
function shockGasModel(gas, variableGamma) {
    if (variableGamma && gas.thermo && !(gas.pinf > 0)) {
        return thermallyPerfectModel(makeThermo(gas.thermo));
    }
    return constantGammaModel(gas);
}

/**
 * Largest shock Mach number a driver can produce (infinite driver pressure):
 * the expansion can only accelerate driver gas to a finite speed.
 */
function maxMachFromModels(drivenModel, driverModel, s1, t4) {
    const uMax = driverModel.maxExpansionVelocity(t4);
    const excess = function (m) { return drivenModel.incidentShock(s1, m).u - uMax; };
    let hi = 2;
    while (excess(hi) < 0 && hi < 1e4) hi *= 2;
    return shockBisect(excess, 1 + 1e-9, hi);
}

function maxMachForDriver(driven, driver, t1, t4, variableGamma) {
    const drivenModel = shockGasModel(driven, variableGamma);
    // The limit does not depend on the driven pressure
    return maxMachFromModels(drivenModel, shockGasModel(driver, variableGamma), drivenModel.state(SHOCK_ATM, t1), t4);
}

/**
 * Compute all shock tube states for a desired incident shock Mach number.
 *
 * @param {Object} input
 * @param {number} input.mach - incident shock Mach number (> 1)
 * @param {number} input.p1 - driven pressure, Pa
 * @param {number} input.t1 - driven temperature, K
 * @param {number} input.t4 - driver temperature, K
 * @param {{gamma:number, mw:number, pinf?:number, thermo?:Array}} input.driven
 * @param {{gamma:number, mw:number, pinf?:number, thermo?:Array}} input.driver
 * @param {number} [input.areaFactor=1] - empirical area-contraction factor
 * @param {function} [input.formatTemperature] - formats a temperature in K
 *   for warning messages (default: kelvin)
 * @param {boolean} [input.variableGamma=false] - use temperature-dependent
 *   specific heats (NASA-9 fits) for gases that have a composition
 * @returns {Object} states 1, 2, 4, 5 (each {p, T, rho, a, u, gamma}) plus
 *   wave speeds and warnings.
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

    const driven = shockGasModel(input.driven, input.variableGamma);
    const driver = shockGasModel(input.driver, input.variableGamma);

    const s1 = driven.state(p1, t1);
    const s2 = driven.incidentShock(s1, m);
    const reflected = driven.reflectedShock(s1, s2);
    const s5 = reflected.state;
    const wr = reflected.speed;

    // State 4: isentropic expansion of the driver from rest to (u2, p2)
    const p4 = driver.expansionSourcePressure(t4, s2.u, s2.p);
    const s4 = driver.state(p4 === null ? p1 : p4, t4);
    const warnings = [];
    let p4psig = null;

    if (p4 === null) {
        s4.p = null;
        s4.rho = null;
        warnings.push({
            code: 'unreachable',
            message: 'Mach ' + m + ' is unreachable with this driver (limit Mach ' +
                maxMachFromModels(driven, driver, s1, t4).toFixed(2) +
                '). Use a lighter or hotter driver gas. Other states are still valid.'
        });
    } else {
        p4psig = (p4 - SHOCK_ATM) / SHOCK_ATM * SHOCK_PSI_PER_ATM / areaFactor;
        if (p4 > EXTREME_DRIVER_PRESSURE_PA && driver.pinf === 0) {
            warnings.push({
                code: 'extreme-driver-pressure',
                message: 'Driver pressure is impractically high (driver limit is Mach ' +
                    maxMachFromModels(driven, driver, s1, t4).toFixed(2) + '). Use a lighter or hotter driver gas.'
            });
        }
    }

    const showT = input.formatTemperature || function (T) { return T.toFixed(0) + ' K'; };
    const hottest = Math.max(s2.T, s5.T);
    if (driven.thermo) {
        if (hottest > Math.min(driven.thermo.tMax, DISSOCIATION_WARNING_K)) {
            warnings.push({
                code: 'high-temperature',
                message: 'Temperatures reach ' + showT(hottest) + '; dissociation and ionization are not modeled.'
            });
        }
    } else if (hottest > HIGH_TEMPERATURE_WARNING_K && driven.pinf === 0) {
        warnings.push({
            code: 'high-temperature',
            message: 'Temperatures reach ' + showT(hottest) + '; constant γ overestimates them above about ' +
                showT(HIGH_TEMPERATURE_WARNING_K) + '. Try the variable-γ gas model.'
        });
    }

    return {
        mach: m,
        gamma: s1.gamma, mw: driven.mw, pinf: driven.pinf,
        variableGamma: !!(driven.thermo || driver.thermo),
        state1: s1,
        state2: s2,
        state4: s4,
        state5: s5,
        shockSpeed: m * s1.a,
        reflectedShockSpeed: wr,
        reflectedMach: (wr + s2.u) / s2.a,
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

    const mMax = maxMachForDriver(input.driven, input.driver, input.t1, input.t4, input.variableGamma);
    // p4(M) increases monotonically from p1 at M = 1 to infinity at mMax
    return shockBisect(function (m) {
        const p4 = computeShockTube({
            mach: m, p1: input.p1, t1: input.t1, t4: input.t4,
            driven: input.driven, driver: input.driver, variableGamma: input.variableGamma
        }).state4.p;
        return p4 === null ? 1 : p4 - input.p4;
    }, 1 + 1e-12, mMax);
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        ShockInputError, computeShockTube, machFromDriverPressure,
        maxMachForDriver, shockGasModel
    };
}
