/**
 * Richtmyer-Meshkov Instability (RMI) Growth Rate Calculator
 * Implements linear and nonlinear models for perturbation growth
 */

/**
 * Calculate linear RMI growth rate (Richtmyer's impulse model)
 * @param {number} eta0 - Initial amplitude (m)
 * @param {number} A - Atwood number
 * @param {number} k - Wavenumber (2π/λ) in rad/m
 * @param {number} deltaU - Velocity jump imparted by shock (m/s)
 * @returns {number} Initial growth rate dη/dt (m/s)
 */
function rmiLinearGrowthRate(eta0, A, k, deltaU) {
    return A * k * deltaU * eta0;
}

/**
 * Calculate linear perturbation amplitude vs time
 * η(t) = η₀·(1 + A·k·Δu·t)
 * @param {number} eta0 - Initial amplitude (m)
 * @param {number} A - Atwood number
 * @param {number} k - Wavenumber (rad/m)
 * @param {number} deltaU - Velocity jump (m/s)
 * @param {number} t - Time after shock (s)
 * @returns {number} Amplitude at time t (m)
 */
function rmiLinearAmplitude(eta0, A, k, deltaU, t) {
    return eta0 * (1.0 + A * k * deltaU * t);
}

/**
 * Calculate nonlinear RMI growth using Zhang-Sohn model
 * Based on Zhang & Sohn (1997) potential flow model
 * Valid for ka₀ << 1 and moderate Atwood numbers
 * @param {number} eta0 - Initial amplitude (m)
 * @param {number} A - Atwood number
 * @param {number} k - Wavenumber (rad/m)
 * @param {number} deltaU - Velocity jump (m/s)
 * @param {number} t - Time after shock (s)
 * @returns {number} Amplitude at time t (m)
 */
function rmiZhangSohnAmplitude(eta0, A, k, deltaU, t) {
    if (Math.abs(A) < 1e-10) {
        return eta0;
    }

    const V0 = rmiLinearGrowthRate(eta0, A, k, deltaU); // Initial growth rate

    // Avoid division issues when V0 is very small
    if (Math.abs(V0) < 1e-12) {
        return eta0;
    }

    const tau = Math.abs(eta0 / V0); // Characteristic time scale
    const tNorm = t / tau;

    // Zhang-Sohn potential flow solution
    // For A > 0 (light pushing heavy): η(t) = η₀ / (1 - t/τ)
    // For A < 0 (heavy pushing light): η(t) = -η₀ / (1 + t/τ)

    // Zhang-Sohn is only valid for t/τ < 0.5 (before nonlinearity saturates)
    // For late times, transition to self-similar turbulent mixing

    if (tNorm < 0.5) {
        // Early-time potential flow regime
        if (A > 0) {
            return eta0 / (1.0 - tNorm);
        } else {
            return -eta0 / (1.0 + tNorm);
        }
    } else if (tNorm < 1.0) {
        // Transition regime (0.5 < t/τ < 1.0)
        const eta_potential = A > 0 ? eta0 / (1.0 - tNorm) : -eta0 / (1.0 + tNorm);

        // Late-time turbulent mixing: h(t) ≈ 0.28·|A|·|Δu|·t
        const h_turb = eta0 + 0.28 * Math.abs(A) * Math.abs(deltaU) * t;

        // Smooth blend factor
        const blend = (tNorm - 0.5) / 0.5; // 0 at t/τ=0.5, 1 at t/τ=1.0
        return eta_potential * (1.0 - blend) + h_turb * blend;
    } else {
        // Late-time turbulent regime (t/τ > 1)
        return eta0 + 0.28 * Math.abs(A) * Math.abs(deltaU) * t;
    }
}

/**
 * Calculate nonlinear RMI growth using Dimonte-Ramaprabhu model
 * Based on buoyancy-drag model from Dimonte & Ramaprabhu (2010)
 * Better for high Atwood numbers and late times
 * @param {number} eta0 - Initial amplitude (m)
 * @param {number} A - Atwood number
 * @param {number} k - Wavenumber (rad/m)
 * @param {number} deltaU - Velocity jump (m/s)
 * @param {number} t - Time after shock (s)
 * @returns {number} Amplitude at time t (m)
 */
function rmiDimonteAmplitude(eta0, A, k, deltaU, t) {
    const V0 = rmiLinearGrowthRate(eta0, A, k, deltaU);

    // Dimonte-Ramaprabhu model parameters
    // h(t) = V₀·t·(1 + C·A·k·V₀·t)^(-1/2)
    // where C ≈ 0.5-1.0 for single-mode RMI (C=2 is too aggressive for drag)
    const C = 0.7;

    if (Math.abs(V0) < 1e-10) {
        return eta0;
    }

    const term = 1 + C * Math.abs(A) * k * Math.abs(V0) * t;
    const h = Math.abs(V0) * t / Math.sqrt(term);

    return eta0 + Math.sign(V0) * h;
}

/**
 * Generate time series of perturbation growth
 * @param {Object} params - {eta0, lambda, A, deltaU, tMax, numPoints, model}
 *   - eta0: initial amplitude (m)
 *   - lambda: wavelength (m)
 *   - A: Atwood number
 *   - deltaU: velocity jump from shock (m/s)
 *   - tMax: maximum time (s)
 *   - numPoints: number of time points
 *   - model: 'linear', 'zhang-sohn', or 'dimonte' (default: 'all')
 * @returns {Object} {times: Array, linear: Array, zhangSohn: Array, dimonte: Array}
 */
function calculateGrowthTimeSeries(params) {
    const {eta0, lambda, A, deltaU, tMax, numPoints = 100, model = 'all'} = params;

    const k = 2 * Math.PI / lambda;
    const times = [];
    const linear = [];
    const zhangSohn = [];
    const dimonte = [];

    for (let i = 0; i < numPoints; i++) {
        const t = (i / (numPoints - 1)) * tMax;
        times.push(t);

        if (model === 'all' || model === 'linear') {
            linear.push(rmiLinearAmplitude(eta0, A, k, deltaU, t));
        }

        if (model === 'all' || model === 'zhang-sohn') {
            zhangSohn.push(rmiZhangSohnAmplitude(eta0, A, k, deltaU, t));
        }

        if (model === 'all' || model === 'dimonte') {
            dimonte.push(rmiDimonteAmplitude(eta0, A, k, deltaU, t));
        }
    }

    return {times, linear, zhangSohn, dimonte};
}

/**
 * Calculate characteristic RMI parameters
 * @param {Object} params - {eta0, lambda, A, deltaU}
 * @returns {Object} Key dimensionless and characteristic parameters
 */
function calculateRMIParameters(params) {
    const {eta0, lambda, A, deltaU} = params;

    const k = 2 * Math.PI / lambda;
    const V0 = rmiLinearGrowthRate(eta0, A, k, deltaU);

    // Initial perturbation parameter
    const ka0 = k * eta0;

    // Characteristic time (time for amplitude to double in linear regime)
    const tDouble = eta0 / Math.abs(V0);

    // Nonlinear time scale (Zhang-Sohn)
    const tNL = 1 / (k * Math.abs(V0));

    return {
        wavenumber: k,
        initialGrowthRate: V0,
        ka0: ka0,
        doublingTime: tDouble,
        nonlinearTimeScale: tNL
    };
}

/**
 * Calculate velocity jump imparted to interface by shock
 * For shock-interface interaction, Δu is the change in interface velocity
 * @param {number} uInterface - Post-shock interface velocity (m/s)
 * @param {number} uInitial - Initial interface velocity (typically 0) (m/s)
 * @returns {number} Velocity jump Δu (m/s)
 */
function calculateVelocityJump(uInterface, uInitial = 0) {
    return uInterface - uInitial;
}

// Export for use in browser and Node.js
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        rmiLinearGrowthRate,
        rmiLinearAmplitude,
        rmiZhangSohnAmplitude,
        rmiDimonteAmplitude,
        calculateGrowthTimeSeries,
        calculateRMIParameters,
        calculateVelocityJump
    };
}
