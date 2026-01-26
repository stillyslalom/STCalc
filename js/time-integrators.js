/**
 * Time Integration Strategies for Euler Solver
 * 
 * This module provides various time integration schemes for solving
 * compressible Euler equations with different stability and accuracy properties.
 */

/**
 * Abstract base class for time integration schemes
 */
class TimeIntegrator {
    /**
     * @param {string} name - Integrator name
     * @param {number} cflRecommended - Recommended CFL number for stability
     */
    constructor(name, cflRecommended) {
        this.name = name;
        this.cflRecommended = cflRecommended;
    }

    /**
     * Take one time step
     * @param {EulerSolver} solver - The Euler solver instance
     * @returns {number} - The time step size used
     */
    step(solver) {
        throw new Error('step() must be implemented by subclass');
    }

    /**
     * Get number of buffer arrays required
     * @returns {number} - Number of Float64Arrays needed
     */
    getRequiredBuffers() {
        throw new Error('getRequiredBuffers() must be implemented by subclass');
    }

    /**
     * Get buffer names for initialization
     * @returns {Array<string>} - Array of buffer names
     */
    getBufferNames() {
        throw new Error('getBufferNames() must be implemented by subclass');
    }

    /**
     * Calculate time step with visualization synchronization
     * Adjusts dt to hit visualization times exactly by reducing timestep
     * when approaching a snapshot time.
     * 
     * @param {EulerSolver} solver - The solver instance
     * @param {number} maxWaveSpeed - Maximum wave speed from flux calculation
     * @returns {number} - The adjusted time step
     */
    calculateTimeStep(solver, maxWaveSpeed) {
        // Normal CFL-based time step
        let dt = solver.cfl * solver.dx / maxWaveSpeed;
        
        // Don't overshoot final time
        if (solver.t + dt > solver.finalTime) {
            dt = solver.finalTime - solver.t;
        }
        
        // Check if we're approaching a visualization time
        const timeToNextViz = solver.nextXTStore - solver.t;
        
        if (timeToNextViz > 0 && timeToNextViz < 2.5 * dt) {
            // We're close to visualization time - adjust dt to hit it exactly
            // Use integer number of substeps for accuracy
            const numSubsteps = Math.max(1, Math.ceil(timeToNextViz / dt));
            dt = timeToNextViz / numSubsteps;
        }
        
        return dt;
    }

    /**
     * Common post-step operations
     * @param {EulerSolver} solver - The Euler solver instance
     */
    postStep(solver) {
        // Update time and step counter
        solver.t += solver.dt;
        solver.timeSteps++;
        
        // Update primitive variables
        solver.updatePrimitives();
        
        // Update Lagrangian tracers
        solver.updateTracers();
        
        // Update gas properties based on interface positions (sharp interface tracking)
        solver.updateGasProperties();
        
        // Store x-t data if we've hit the target time (within numerical tolerance)
        if (Math.abs(solver.t - solver.nextXTStore) < 1e-10 || solver.t > solver.nextXTStore) {
            solver.storeXTData();
            solver.nextXTStore += solver.xtInterval;
        }
    }
}

/**
 * Second-order Runge-Kutta (RK2) Time Integrator
 * 
 * Also known as the midpoint method or Heun's method variant.
 * This is a two-stage, second-order accurate explicit method.
 * 
 * Scheme:
 *   U* = U^n + dt * L(U^n)
 *   U^{n+1} = 0.5 * U^n + 0.5 * (U* + dt * L(U*))
 * 
 * Properties:
 * - Order of accuracy: 2
 * - Stages: 2
 * - Recommended CFL: 0.4
 * - Memory: 3 buffers (U, Un, Uk)
 */
class RK2Integrator extends TimeIntegrator {
    constructor() {
        super('RK2', 0.4);
    }

    getRequiredBuffers() {
        return 3;
    }

    getBufferNames() {
        return ['U', 'Un', 'Uk'];
    }

    step(solver) {
        // Stage 1: compute flux and update
        const maxWaveSpeed = solver.computeFluxes();
        
        // Calculate time step with visualization synchronization
        solver.dt = this.calculateTimeStep(solver, maxWaveSpeed);
        
        // RK Stage 1: U* = U^n + dt * L(U^n)
        solver.updateConservative(solver.U, solver.Uk, solver.dt);
        
        // Stage 2: compute flux from U*
        // Copy Uk to U temporarily for flux computation
        const n = solver.nx * 3;
        for (let i = 0; i < n; i++) {
            solver.Un[i] = solver.U[i];
            solver.U[i] = solver.Uk[i];
        }
        
        solver.computeFluxes();
        
        // RK Stage 2: U^{n+1} = 0.5 * U^n + 0.5 * (U* + dt * L(U*))
        solver.updateConservative(solver.Uk, solver.Uk, solver.dt);
        
        for (let i = 0; i < n; i++) {
            solver.U[i] = 0.5 * solver.Un[i] + 0.5 * solver.Uk[i];
        }
        
        // Post-step operations
        this.postStep(solver);
        
        return solver.dt;
    }
}

/**
 * Strong Stability Preserving (SSP) Runge-Kutta Integrator
 * 
 * Implements the optimal 4-stage, 3rd-order SSP-RK scheme SSPRK(4,3).
 * This scheme maintains the total variation diminishing (TVD) property
 * and provides enhanced stability for hyperbolic conservation laws.
 * 
 * Reference:
 *   Gottlieb, S., Shu, C.-W., & Tadmor, E. (2001).
 *   "Strong Stability-Preserving High-Order Time Discretization Methods"
 *   SIAM Review, 43(1), 89-112.
 * 
 * Scheme (in Shu-Osher form):
 *   U^(1) = U^n + 0.5 * dt * L(U^n)
 *   U^(2) = U^(1) + 0.5 * dt * L(U^(1))
 *   U^(3) = (2/3)*U^n + (1/6)*U^(2) + (1/6)*dt*L(U^(2))
 *   U^{n+1} = 0.5*U^(3) + 0.5*U^(3) + 0.5*dt*L(U^(3))
 * 
 * Properties:
 * - Order of accuracy: 3
 * - Stages: 4
 * - Recommended CFL: 0.8 (2x larger than RK2)
 * - Memory: 3 buffers (U, Un, Uk) - optimized for reuse
 * - SSP coefficient: c = 2.0
 * 
 * Note: While 4 stages are used, the scheme is designed to be SSP with
 * CFL up to 2.0 times the forward Euler limit, making it more efficient
 * than RK2 for many problems despite the extra stages.
 */
class SSPIntegrator extends TimeIntegrator {
    constructor() {
        super('SSP', 0.8);
    }

    getRequiredBuffers() {
        return 3;  // Optimized to use same as RK2
    }

    getBufferNames() {
        return ['U', 'Un', 'Uk'];
    }

    step(solver) {
        // Compute initial flux and time step
        const maxWaveSpeed = solver.computeFluxes();
        
        // Calculate time step with visualization synchronization
        solver.dt = this.calculateTimeStep(solver, maxWaveSpeed);
        
        const n = solver.nx * 3;
        
        // Save initial state in Un
        for (let i = 0; i < n; i++) {
            solver.Un[i] = solver.U[i];
        }
        
        // Stage 1: U^(1) = U^n + 0.5 * dt * L(U^n)
        solver.updateConservative(solver.U, solver.Uk, solver.dt);
        for (let i = 0; i < n; i++) {
            solver.U[i] = solver.U[i] + 0.5 * (solver.Uk[i] - solver.U[i]);
        }
        
        // Stage 2: U^(2) = U^(1) + 0.5 * dt * L(U^(1))
        solver.computeFluxes();
        solver.updateConservative(solver.U, solver.Uk, solver.dt);
        for (let i = 0; i < n; i++) {
            solver.U[i] = solver.U[i] + 0.5 * (solver.Uk[i] - solver.U[i]);
        }
        
        // Stage 3: U^(3) = (2/3)*U^n + (1/6)*U^(2) + (1/6)*dt*L(U^(2))
        solver.computeFluxes();
        solver.updateConservative(solver.U, solver.Uk, solver.dt);
        for (let i = 0; i < n; i++) {
            solver.U[i] = (2.0/3.0) * solver.Un[i] + (1.0/6.0) * solver.U[i] + (1.0/6.0) * solver.Uk[i];
        }
        
        // Stage 4: U^{n+1} = 0.5*U^(3) + 0.5*(U^(3) + dt*L(U^(3)))
        solver.computeFluxes();
        solver.updateConservative(solver.U, solver.Uk, solver.dt);
        for (let i = 0; i < n; i++) {
            solver.U[i] = 0.5 * solver.U[i] + 0.5 * solver.Uk[i];
        }
        
        // Post-step operations
        this.postStep(solver);
        
        return solver.dt;
    }
}

/**
 * Factory for creating time integrators
 */
class IntegratorFactory {
    static create(name) {
        switch (name) {
            case 'RK2':
                return new RK2Integrator();
            case 'SSP':
                return new SSPIntegrator();
            default:
                console.warn(`Unknown integrator: ${name}, defaulting to RK2`);
                return new RK2Integrator();
        }
    }

    static getAvailable() {
        return [
            { name: 'RK2', cfl: 0.4, description: '2nd-order Runge-Kutta' },
            { name: 'SSP', cfl: 0.8, description: '3rd-order Strong Stability Preserving' }
        ];
    }
}
