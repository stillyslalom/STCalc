/**
 * Exact Riemann Solver for compressible gas dynamics
 * Solves for interface conditions between two ideal or stiffened gases with different properties
 */

/**
 * Solve the exact Riemann problem for two stiffened gases.
 *
 * Each side obeys p + pinf = rho * R * T with a^2 = gamma * (p + pinf) / rho.
 * Setting pinf = 0 (or omitting it) gives an ideal gas. All ideal-gas wave
 * relations hold in terms of the shifted pressure p + pinf.
 *
 * @param {Object} leftState - {p: pressure (Pa), T: temperature (K), u: velocity (m/s), gamma, mw (g/mol), pinf (Pa, optional)}
 * @param {Object} rightState - {p, T, u, gamma, mw, pinf}
 * @param {number} ru - Universal gas constant (J/kmol-K), default 8314.51
 * @returns {Object} {pStar, uStar, leftPost, rightPost, waveConfig, converged}
 *   - pStar: interface pressure (Pa)
 *   - uStar: interface velocity (m/s)
 *   - leftPost: {p, T, rho, u, a} post-wave state on left side
 *   - rightPost: {p, T, rho, u, a} post-wave state on right side
 *   - waveConfig: {left: 'shock'|'rarefaction', right: 'shock'|'rarefaction'}
 */
function solveExactRiemann(leftState, rightState, ru = 8314.51) {
    const side = function (state) {
        const R = ru / state.mw;
        const pinf = state.pinf || 0;
        const pt = state.p + pinf;          // shifted pressure
        return {
            R, pinf, pt,
            p: state.p,
            T: state.T,
            u: state.u || 0,
            gamma: state.gamma,
            rho: pt / (R * state.T),
            a: Math.sqrt(state.gamma * R * state.T)
        };
    };
    const L = side(leftState);
    const R = side(rightState);

    // Lowest admissible pressure: both shifted pressures must stay positive
    const pFloor = -Math.min(L.pinf, R.pinf);
    const scale = Math.max(L.pt, R.pt);

    // f(p) = fL + fR + (uR - uL) is monotonically increasing and concave
    const residual = function (p) {
        return pressureFunction(p + L.pinf, L.pt, L.rho, L.a, L.gamma) +
               pressureFunction(p + R.pinf, R.pt, R.rho, R.a, R.gamma) +
               (R.u - L.u);
    };
    const slope = function (p) {
        return pressureFunctionDerivative(p + L.pinf, L.pt, L.rho, L.a, L.gamma) +
               pressureFunctionDerivative(p + R.pinf, R.pt, R.rho, R.a, R.gamma);
    };

    // Acoustic (linearised) initial guess, valid for any equation of state
    const zL = L.rho * L.a, zR = R.rho * R.a;
    let pStar = (zR * L.p + zL * R.p - zL * zR * (R.u - L.u)) / (zL + zR);
    const pMin = pFloor + 1e-12 * scale;
    if (!isFinite(pStar) || pStar <= pMin) pStar = pMin;

    const maxIter = 200;
    const tol = 1e-9 * (L.a + R.a);
    let converged = false;

    for (let iter = 0; iter < maxIter; iter++) {
        const f = residual(pStar);
        if (Math.abs(f) < tol) {
            converged = true;
            break;
        }
        let next = pStar - f / slope(pStar);
        // Keep the iterate admissible; halve towards the floor if Newton overshoots
        if (!(next > pMin)) next = 0.5 * (pStar + pMin);
        if (Math.abs(next - pStar) <= 1e-14 * scale) {
            pStar = next;
            converged = true;
            break;
        }
        pStar = next;
    }

    if (!converged) {
        console.warn('Riemann solver did not converge, using best estimate');
    }

    const fL = pressureFunction(pStar + L.pinf, L.pt, L.rho, L.a, L.gamma);
    const fR = pressureFunction(pStar + R.pinf, R.pt, R.rho, R.a, R.gamma);
    const uStar = 0.5 * (L.u + R.u) + 0.5 * (fR - fL);

    const waveConfig = {
        left: pStar > L.p ? 'shock' : 'rarefaction',
        right: pStar > R.p ? 'shock' : 'rarefaction'
    };

    return {
        pStar,
        uStar,
        leftPost: calculatePostWaveState(L, pStar, uStar, waveConfig.left),
        rightPost: calculatePostWaveState(R, pStar, uStar, waveConfig.right),
        waveConfig,
        converged
    };
}

/**
 * Pressure function f(p) for Riemann solver: velocity change across the wave.
 * p and p0 are shifted pressures (p + pinf).
 * @private
 */
function pressureFunction(p, p0, rho0, a0, gamma) {
    if (p > p0) {
        // Shock wave
        const A = 2 / ((gamma + 1) * rho0);
        const B = (gamma - 1) / (gamma + 1) * p0;
        return (p - p0) * Math.sqrt(A / (p + B));
    } else {
        // Rarefaction wave
        const exponent = (gamma - 1) / (2 * gamma);
        return (2 * a0 / (gamma - 1)) * (Math.pow(p / p0, exponent) - 1);
    }
}

/**
 * Derivative of pressure function f'(p)
 * @private
 */
function pressureFunctionDerivative(p, p0, rho0, a0, gamma) {
    if (p > p0) {
        // Shock wave
        const A = 2 / ((gamma + 1) * rho0);
        const B = (gamma - 1) / (gamma + 1) * p0;
        const sqrtTerm = Math.sqrt(A / (p + B));
        return sqrtTerm * (1 - (p - p0) / (2 * (p + B)));
    } else {
        // Rarefaction wave
        const exponent = (gamma - 1) / (2 * gamma);
        return (1 / (rho0 * a0)) * Math.pow(p / p0, -((gamma + 1) / (2 * gamma)));
    }
}

/**
 * Calculate post-wave state (after shock or rarefaction)
 * @private
 */
function calculatePostWaveState(initialState, pStar, uStar, waveType) {
    const {pt: pt0, rho: rho0, gamma, R, pinf} = initialState;
    const ratio = (pStar + pinf) / pt0;     // shifted pressure ratio
    let rho;

    if (waveType === 'shock') {
        const gp = gamma + 1;
        const gm = gamma - 1;
        rho = rho0 * (gp * ratio + gm) / (gm * ratio + gp);
    } else {
        rho = rho0 * Math.pow(ratio, 1 / gamma);
    }

    const T = (pStar + pinf) / (rho * R);
    return {
        p: pStar,
        rho,
        T,
        u: uStar,
        a: Math.sqrt(gamma * R * T)
    };
}

/**
 * Calculate Atwood number for two-gas interface
 * @param {number} rho1 - Density of gas 1 (kg/m³)
 * @param {number} rho2 - Density of gas 2 (kg/m³)
 * @returns {number} Atwood number A = (rho2 - rho1) / (rho2 + rho1)
 */
function calculateAtwoodNumber(rho1, rho2) {
    return (rho2 - rho1) / (rho2 + rho1);
}

/**
 * Solve shocked-interface problem for RMI experiments
 * Scenario: shock propagates through light gas and impacts light/heavy interface
 * @param {Object} shockedLightState - State of light gas after shock passage {p, T, u, gamma, mw, pinf}
 * @param {Object} ambientHeavyState - Ambient state of heavy gas before shock arrival {p, T, u, gamma, mw, pinf}
 * @param {number} ru - Universal gas constant, default 8314.51
 * @returns {Object} Interface conditions and Atwood number
 */
function solveShockedInterface(shockedLightState, ambientHeavyState, ru = 8314.51) {
    const solution = solveExactRiemann(shockedLightState, ambientHeavyState, ru);

    // Calculate Atwood number (using post-shock densities at interface)
    const A = calculateAtwoodNumber(solution.leftPost.rho, solution.rightPost.rho);

    return {
        interfaceVelocity: solution.uStar,
        interfacePressure: solution.pStar,
        lightGasPost: solution.leftPost,
        heavyGasPost: solution.rightPost,
        atwoodNumber: A,
        waveConfig: solution.waveConfig,
        converged: solution.converged
    };
}

// Export for use in browser and Node.js
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        solveExactRiemann,
        calculateAtwoodNumber,
        solveShockedInterface
    };
}
