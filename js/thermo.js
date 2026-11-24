/**
 * Thermally perfect gas properties: an ideal gas (p = rho R T) whose specific
 * heats depend on temperature, built from NASA-9 fits (js/nasa9-data.js).
 * Composition is frozen: no dissociation, ionization or reaction.
 *
 * Pure functions only: no DOM access.
 */

const THERMO_RU = 8314.51;   // J/(kmol K), as in the rest of the calculator

/**
 * Molar properties of one species from its NASA-9 fit.
 * Outside the fitted range the heat capacity is held at its end value, which
 * keeps enthalpy and entropy smooth and monotonic (e.g. for a driver gas that
 * expands to very low temperature).
 * @returns {{mw:number, tMin:number, tMax:number, cpR:function, hR:function, sR:function}}
 *   cpR = cp/R, hR = h/R (in K), sR = s°/R, all as functions of T (K)
 */
function nasa9Species(name) {
    const record = NASA9_DATA[name];
    if (!record) throw new Error('No NASA-9 data for species ' + name);
    const intervals = record.intervals;
    const tMin = intervals[0].Tmin;
    const tMax = intervals[intervals.length - 1].Tmax;

    const pick = function (T) {
        for (let i = 0; i < intervals.length; i++) {
            if (T <= intervals[i].Tmax) return intervals[i];
        }
        return intervals[intervals.length - 1];
    };
    const cpFit = function (T) {
        const a = pick(T).a;
        return a[0] / (T * T) + a[1] / T + a[2] + T * (a[3] + T * (a[4] + T * (a[5] + T * a[6])));
    };
    const hFit = function (T) {
        const iv = pick(T), a = iv.a;
        return -a[0] / T + a[1] * Math.log(T) + iv.b1 +
               T * (a[2] + T * (a[3] / 2 + T * (a[4] / 3 + T * (a[5] / 4 + T * a[6] / 5))));
    };
    const sFit = function (T) {
        const iv = pick(T), a = iv.a;
        return -a[0] / (2 * T * T) - a[1] / T + a[2] * Math.log(T) + iv.b2 +
               T * (a[3] + T * (a[4] / 2 + T * (a[5] / 3 + T * a[6] / 4)));
    };
    const ends = [tMin, tMax].map(function (T) { return { T: T, cp: cpFit(T), h: hFit(T), s: sFit(T) }; });
    const outside = function (T) { return T < tMin ? ends[0] : (T > tMax ? ends[1] : null); };

    return {
        mw: record.mw, tMin: tMin, tMax: tMax,
        cpR: function (T) { const e = outside(T); return e ? e.cp : cpFit(T); },
        hR: function (T) { const e = outside(T); return e ? e.h + e.cp * (T - e.T) : hFit(T); },
        sR: function (T) { const e = outside(T); return e ? e.s + e.cp * Math.log(T / e.T) : sFit(T); }
    };
}

/** A species with constant specific heats (used for user-defined gases). */
function constantCpSpecies(gamma, mw) {
    const cp = gamma / (gamma - 1);
    return {
        mw: mw, tMin: 0, tMax: Infinity,
        cpR: function () { return cp; },
        hR: function (T) { return cp * T; },
        sR: function (T) { return cp * Math.log(T); }
    };
}

/**
 * Build the property functions of a frozen mixture.
 * @param {Array} composition - entries {species: 'N2', x} or {gamma, mw, x},
 *   x being mole fractions (normalized here)
 * @returns {Object} per-unit-mass functions of temperature:
 *   R, mw, cp(T), h(T), s0(T) (entropy at unit pressure), gamma(T), a(T),
 *   temperatureFromEnthalpy(h, guess), and the fitted range tMin..tMax
 */
function makeThermo(composition) {
    const parts = composition.map(function (c) {
        return { x: c.x, sp: c.species ? nasa9Species(c.species) : constantCpSpecies(c.gamma, c.mw) };
    });
    const total = parts.reduce(function (sum, p) { return sum + p.x; }, 0);
    parts.forEach(function (p) { p.x /= total; });

    const mw = parts.reduce(function (sum, p) { return sum + p.x * p.sp.mw; }, 0);
    const R = THERMO_RU / mw;
    const molar = function (fn, T) {
        let sum = 0;
        for (let i = 0; i < parts.length; i++) sum += parts[i].x * parts[i].sp[fn](T);
        return sum;
    };
    const cp = function (T) { return R * molar('cpR', T); };
    const h = function (T) { return R * molar('hR', T); };
    const gamma = function (T) { const c = molar('cpR', T); return c / (c - 1); };

    return {
        mw: mw, R: R, cp: cp, h: h, gamma: gamma,
        s0: function (T) { return R * molar('sR', T); },
        a: function (T) { return Math.sqrt(gamma(T) * R * T); },
        tMin: Math.max.apply(null, parts.map(function (p) { return p.sp.tMin; })),
        tMax: Math.min.apply(null, parts.map(function (p) { return p.sp.tMax; })),
        /** Invert h(T); h is strictly increasing because cp > 0. */
        temperatureFromEnthalpy: function (target, guess) {
            let T = guess > 0 ? guess : 300;
            for (let i = 0; i < 100; i++) {
                const step = (target - h(T)) / cp(T);
                const next = Math.max(0.2 * T, T + step);   // never jump to T <= 0
                if (Math.abs(next - T) <= 1e-12 * T) return next;
                T = next;
            }
            return T;
        }
    };
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { nasa9Species, constantCpSpecies, makeThermo };
}
