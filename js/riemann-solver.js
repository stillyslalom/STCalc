/**
 * Exact Riemann Solver for compressible gas dynamics
 * Solves for interface conditions between two ideal gases with different properties
 */

/**
 * Solve the exact Riemann problem for two ideal gases
 * @param {Object} leftState - {p: pressure (Pa), T: temperature (K), u: velocity (m/s), gamma: ratio of specific heats, mw: molecular weight (g/mol)}
 * @param {Object} rightState - {p, T, u, gamma, mw}
 * @param {number} ru - Universal gas constant (J/mol-K), default 8314.51
 * @returns {Object} {pStar, uStar, leftPost, rightPost, waveConfig}
 *   - pStar: interface pressure (Pa)
 *   - uStar: interface velocity (m/s)
 *   - leftPost: {p, T, rho, u, a} post-wave state on left side
 *   - rightPost: {p, T, rho, u, a} post-wave state on right side
 *   - waveConfig: {left: 'shock'|'rarefaction', right: 'shock'|'rarefaction'}
 */
function solveExactRiemann(leftState, rightState, ru = 8314.51) {
    // Calculate gas constants
    const RL = ru / leftState.mw;
    const RR = ru / rightState.mw;

    // Initial velocities (default to 0 if not provided)
    const uL0 = leftState.u || 0;
    const uR0 = rightState.u || 0;

    // Initial densities
    const rhoL = leftState.p / (RL * leftState.T);
    const rhoR = rightState.p / (RR * rightState.T);

    // Sound speeds
    const aL = Math.sqrt(leftState.gamma * RL * leftState.T);
    const aR = Math.sqrt(rightState.gamma * RR * rightState.T);

    // Pressure ratio guess (use acoustic approximation)
    const pRatio = leftState.p / rightState.p;
    let pStar;

    if (pRatio > 1) {
        // Shock likely on right
        pStar = Math.max(leftState.p, rightState.p);
    } else {
        // Shock likely on left
        pStar = Math.min(leftState.p, rightState.p);
    }

    // Better initial guess using two-rarefaction approximation
    const gamL = leftState.gamma;
    const gamR = rightState.gamma;

    // Account for initial velocity difference in guess
    const duInit = uR0 - uL0;
    const pGuess = Math.pow(
        (aL + aR - 0.5 * (gamL - 1) * duInit) /
        (aL / Math.pow(leftState.p, (gamL - 1) / (2 * gamL)) + aR / Math.pow(rightState.p, (gamR - 1) / (2 * gamR))),
        (2 * gamL * gamR) / (gamL - 1 + gamR - 1)
    );

    if (!isNaN(pGuess) && pGuess > 0) {
        pStar = pGuess;
    }

    // Newton-Raphson iteration to find p*
    const maxIter = 100;
    const tol = 1e-6;
    let iter = 0;

    while (iter < maxIter) {
        // Calculate f(p*) and f'(p*) for both sides
        const fL = pressureFunction(pStar, leftState.p, rhoL, aL, gamL);
        const fR = pressureFunction(pStar, rightState.p, rhoR, aR, gamR);
        const dfL = pressureFunctionDerivative(pStar, leftState.p, rhoL, aL, gamL);
        const dfR = pressureFunctionDerivative(pStar, rightState.p, rhoR, aR, gamR);

        // The condition is: f = fL + fR + (uR - uL) = 0
        const f = fL + fR + (uR0 - uL0);
        const df = dfL + dfR;

        if (Math.abs(f) < tol) {
            break;
        }

        // Newton-Raphson update
        const dp = -f / df;
        pStar = pStar + dp;

        // Ensure pressure stays positive
        pStar = Math.max(pStar, 1e-6);

        iter++;
    }

    if (iter >= maxIter) {
        console.warn('Riemann solver did not converge, using best estimate');
    }

    // Calculate u* from p* (velocity in star region)
    // From the Julia code: u_star = 0.5 * (u_L + u_R) + 0.5 * (f_R - f_L)
    const fL = velocityFunction(pStar, leftState.p, rhoL, aL, gamL);
    const fR = velocityFunction(pStar, rightState.p, rhoR, aR, gamR);
    const uStar = 0.5 * (uL0 + uR0) + 0.5 * (fR - fL);

    // Determine wave configuration
    const waveConfig = {
        left: pStar > leftState.p ? 'shock' : 'rarefaction',
        right: pStar > rightState.p ? 'shock' : 'rarefaction'
    };

    // Calculate post-wave states
    const leftPost = calculatePostWaveState(
        {p: leftState.p, rho: rhoL, u: uL0, a: aL, gamma: gamL, R: RL, T: leftState.T},
        pStar,
        uStar,
        waveConfig.left
    );

    const rightPost = calculatePostWaveState(
        {p: rightState.p, rho: rhoR, u: uR0, a: aR, gamma: gamR, R: RR, T: rightState.T},
        pStar,
        uStar,
        waveConfig.right
    );

    return {
        pStar,
        uStar,
        leftPost,
        rightPost,
        waveConfig
    };
}

/**
 * Pressure function f(p) for Riemann solver
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
 * Calculate particle velocity as function of pressure
 * @private
 */
function velocityFunction(p, p0, rho0, a0, gamma) {
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
 * Calculate post-wave state (after shock or rarefaction)
 * @private
 */
function calculatePostWaveState(initialState, pStar, uStar, waveType) {
    const {p: p0, rho: rho0, a: a0, gamma, R, T: T0} = initialState;

    let rho, T, a;

    if (waveType === 'shock') {
        // Shock relations
        const pRatio = pStar / p0;
        const gp = gamma + 1;
        const gm = gamma - 1;

        rho = rho0 * (gp * pRatio + gm) / (gm * pRatio + gp);
        T = pStar / (rho * R);
        a = Math.sqrt(gamma * R * T);
    } else {
        // Rarefaction relations (isentropic)
        const pRatio = pStar / p0;
        const exponent = (gamma - 1) / gamma;

        rho = rho0 * Math.pow(pRatio, 1 / gamma);
        T = T0 * Math.pow(pRatio, exponent);
        a = Math.sqrt(gamma * R * T);
    }

    return {
        p: pStar,
        rho,
        T,
        u: uStar,
        a
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
 * @param {Object} shockedLightState - State of light gas after shock passage {p, T, gamma, mw}
 * @param {Object} ambientHeavyState - Ambient state of heavy gas before shock arrival {p, T, gamma, mw}
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
        waveConfig: solution.waveConfig
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
