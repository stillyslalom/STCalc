/**
 * Gas database and equation-of-state helpers.
 *
 * Every material is described by the stiffened-gas equation of state
 *
 *     p = (gamma - 1) * rho * e - gamma * pinf
 *     p + pinf = rho * R * T,          R = Ru / mw
 *     a^2 = gamma * (p + pinf) / rho = gamma * R * T
 *
 * With pinf = 0 this is the calorically perfect ideal gas. With pinf > 0 it
 * models a stiff liquid such as water; "mw" is then an effective value chosen
 * so that the reference density is reproduced, not a true molecular weight.
 *
 * This file is the single source of truth for gas properties. It has no DOM
 * dependencies so that it can be unit tested under Node.
 */

const RU = 8314.51;             // Universal gas constant, J/(kmol K)
const STP_REF_T = 293.15;       // Reference state used to define liquids, K
const STP_REF_P = 101325;       // Pa

/** Effective molecular weight that reproduces rhoRef at (pRef, tRef). */
function stiffenedEffectiveMw(pinf, rhoRef, tRef = STP_REF_T, pRef = STP_REF_P) {
    return RU * rhoRef * tRef / (pRef + pinf);
}

// gamma and mw are the constant-gamma values of the original calculator.
// "nasa" lists [species, mole fraction] pairs in js/nasa9-data.js, used by the
// temperature-dependent gas model.
const GAS_DATABASE = {
    "1": { name: "Air", gamma: 1.4020, mw: 28.9660,
           nasa: [["N2", 0.78084], ["O2", 0.209476], ["Ar", 0.009365], ["CO2", 0.000319]] },
    "2": { name: "Argon", gamma: 1.6700, mw: 39.948, nasa: [["Ar", 1]] },
    "3": { name: "Carbon Dioxide", gamma: 1.297, mw: 44.01, nasa: [["CO2", 1]] },
    "4": { name: "Helium", gamma: 1.667, mw: 4.002602, nasa: [["He", 1]] },
    "5": { name: "Neon", gamma: 1.667, mw: 20.1797, nasa: [["Ne", 1]] },
    "6": { name: "Nitrogen", gamma: 1.401, mw: 28.0134, nasa: [["N2", 1]] },
    "7": { name: "Xenon", gamma: 1.667, mw: 131.29, nasa: [["Xe", 1]] },
    "8": { name: "Sulfur Hexafluoride", gamma: 1.092, mw: 146.06, nasa: [["SF6", 1]] },
    "9": { name: "Helium-Neon", gamma: 1.668, mw: 21.98 },
    "10": { name: "Acetone Vapor", gamma: 1.14, mw: 58.08, nasa: [["C3H6O,acetone", 1]] },
    "11": { name: "Hydrogen", gamma: 1.405, mw: 2.01588, nasa: [["H2", 1]] },
    "12": { name: "Krypton", gamma: 1.667, mw: 83.798, nasa: [["Kr", 1]] },
    "13": { name: "Oxygen", gamma: 1.395, mw: 31.9988, nasa: [["O2", 1]] },
    "14": { name: "Methane", gamma: 1.304, mw: 16.0425, nasa: [["CH4", 1]] },
    // Stiffened-gas fit for liquid water: a = 1450 m/s and rho = 998.2 kg/m3
    // at 293.15 K, 1 atm.
    "water": {
        name: "Water (liquid, stiffened gas)",
        gamma: 6.12,
        pinf: 3.43e8,
        mw: stiffenedEffectiveMw(3.43e8, 998.2)
    }
};

/** Normalize a gas record so that pinf is always a finite number. */
function normalizeGas(gas) {
    return {
        name: gas.name,
        gamma: gas.gamma,
        mw: gas.mw,
        pinf: gas.pinf > 0 ? gas.pinf : 0,
        nasa: gas.nasa || null
    };
}

/**
 * Composition of a gas for the temperature-dependent model (see thermo.js):
 * NASA-9 species where known, otherwise one constant-cp pseudo-species.
 * Liquids have no such model and return null.
 */
function gasComposition(gas) {
    if (isStiffened(gas)) return null;
    if (gas.nasa) return gas.nasa.map(function (entry) { return { species: entry[0], x: entry[1] }; });
    return [{ gamma: gas.gamma, mw: gas.mw, x: 1 }];
}

function isStiffened(gas) {
    return !!gas && gas.pinf > 0;
}

function gasConstant(gas) {
    return RU / gas.mw;
}

function eosDensity(gas, p, T) {
    return (p + (gas.pinf || 0)) / (gasConstant(gas) * T);
}

function eosSoundSpeed(gas, T) {
    return Math.sqrt(gas.gamma * gasConstant(gas) * T);
}

/** Validate user-supplied gas properties. Returns an error string or null. */
function validateGasDefinition(def) {
    if (!def.name || !String(def.name).trim()) return 'Please enter a name.';
    if (!isFinite(def.gamma) || def.gamma <= 1) return 'Gamma must be a number greater than 1.';
    if (def.pinf !== undefined && (!isFinite(def.pinf) || def.pinf < 0)) return 'p∞ must be zero or positive.';
    if (!isFinite(def.mw) || def.mw <= 0) {
        return def.pinf > 0 ? 'Reference density must be a positive number.'
                            : 'Molecular weight must be a positive number.';
    }
    return null;
}

/**
 * Properties of a mole-fraction mixture of ideal gases.
 * @param {Array<{gasId:string, fraction:number}>} components
 * @param {Object} gases - lookup of gasId -> gas record
 * @returns {{gamma:number, mw:number, pinf:number, name:string, thermo:Array|null}}
 * @throws {Error} with a user-readable message when the mixture is invalid
 */
function mixtureProperties(components, gases) {
    const valid = [];
    let total = 0;
    for (const c of components || []) {
        const fraction = parseFloat(c.fraction);
        if (!c.gasId || !(fraction > 0)) continue;
        const gas = gases[c.gasId];
        if (!gas) throw new Error('Mixture contains an unknown gas.');
        valid.push({ gas, fraction });
        total += fraction;
    }
    if (valid.length === 0) throw new Error('No gas components specified in the mixture.');
    if (Math.abs(total - 1) > 0.001) {
        throw new Error('Mole fractions must sum to 1.0 (100%). Current total: ' + (total * 100).toFixed(2) + '%.');
    }
    if (valid.length > 1 && valid.some(v => isStiffened(v.gas))) {
        throw new Error('Stiffened-gas liquids cannot be blended in a mole-fraction mixture. Use the liquid on its own.');
    }

    let mw = 0, cv = 0, cp = 0;
    for (const { gas, fraction } of valid) {
        mw += fraction * gas.mw;
        // Molar specific heats in units of Ru
        const cvI = 1 / (gas.gamma - 1);
        cv += fraction * cvI;
        cp += fraction * gas.gamma * cvI;
    }
    // Composition for the temperature-dependent model
    let thermo = [];
    for (const { gas, fraction } of valid) {
        const parts = gasComposition(gas);
        if (!parts) {
            thermo = null;
            break;
        }
        for (const part of parts) thermo.push({ ...part, x: part.x * fraction });
    }

    return {
        thermo: thermo,
        gamma: cp / cv,
        mw: mw / total,
        pinf: valid.length === 1 ? (valid[0].gas.pinf || 0) : 0,
        name: valid.map(v => Math.round(v.fraction * 100) + '% ' + v.gas.name).join(' + ')
    };
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        RU, STP_REF_T, STP_REF_P, GAS_DATABASE,
        stiffenedEffectiveMw, normalizeGas, isStiffened, gasConstant, gasComposition,
        eosDensity, eosSoundSpeed, validateGasDefinition, mixtureProperties
    };
}
