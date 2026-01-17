/**
 * 1D Compressible Euler Solver with HLLC Flux
 * For multicomponent shock tube simulations
 */

class EulerSolver {
    constructor(config) {
        this.nx = config.nx || 500;  // Number of grid cells
        this.cfl = config.cfl || 0.4;  // CFL number for stability
        this.finalTime = config.finalTime || 0.02;  // Final simulation time (seconds)
        
        // Initialize arrays using typed arrays for performance
        this.x = new Float64Array(this.nx);  // Cell centers
        this.dx = 0;  // Grid spacing (set during initialization)
        
        // Conservative variables: [rho, rho*u, E]
        this.U = new Float64Array(this.nx * 3);
        this.Un = new Float64Array(this.nx * 3);  // Next time step
        this.Uk = new Float64Array(this.nx * 3);  // RK stage
        
        // Primitive variables for output
        this.rho = new Float64Array(this.nx);
        this.u = new Float64Array(this.nx);
        this.p = new Float64Array(this.nx);
        this.T = new Float64Array(this.nx);
        
        // Gas properties per cell
        this.gamma = new Float64Array(this.nx);
        this.mw = new Float64Array(this.nx);  // Molecular weight
        this.gasId = new Array(this.nx);  // Gas identifier for each cell
        
        // Flux array
        this.flux = new Float64Array((this.nx + 1) * 3);
        
        // Time tracking
        this.t = 0;
        this.dt = 0;
        this.timeSteps = 0;
        
        // X-T diagram storage
        this.xtData = [];  // Will store {t, x, p} at specified intervals
        this.xtInterval = config.xtInterval || 0.0001;  // Store every 0.1 ms
        this.nextXTStore = 0;
        
        // Lagrangian tracers
        this.tracers = [];
        
        // Regions for sharp interface tracking
        this.regions = [];  // Store {gamma, mw, gasId} for each region between interfaces
        
        // Constants
        this.Ru = 8314.51;  // Universal gas constant J/(kmol·K)
    }
    
    /**
     * Initialize the shock tube with gas slabs
     * @param {Array} slabs - Array of {gas, pressure, temperature, length, position}
     */
    initialize(slabs) {
        // Calculate total length
        const totalLength = slabs.reduce((sum, slab) => sum + slab.length, 0);
        this.dx = totalLength / this.nx;
        
        // Create grid points (cell centers)
        for (let i = 0; i < this.nx; i++) {
            this.x[i] = (i + 0.5) * this.dx;
        }
        
        // Assign gas properties and initial conditions to each cell
        let currentPos = 0;
        let slabIndex = 0;
        
        for (let i = 0; i < this.nx; i++) {
            const cellCenter = this.x[i];
            
            // Find which slab this cell belongs to
            while (slabIndex < slabs.length - 1 && cellCenter >= currentPos + slabs[slabIndex].length) {
                currentPos += slabs[slabIndex].length;
                slabIndex++;
            }
            
            const slab = slabs[slabIndex];
            
            // Store gas properties
            this.gamma[i] = slab.gamma;
            this.mw[i] = slab.mw;
            this.gasId[i] = slab.gasId;
            
            // Calculate gas constant for this gas
            const R = this.Ru / slab.mw;
            
            // Set primitive variables
            const rho = slab.pressure / (R * slab.temperature);
            const u = 0;  // Initially at rest
            const p = slab.pressure;
            
            // Convert to conservative variables
            const E = p / (slab.gamma - 1) + 0.5 * rho * u * u;
            
            const idx = i * 3;
            this.U[idx] = rho;
            this.U[idx + 1] = rho * u;
            this.U[idx + 2] = E;
        }
        
        // Store region properties for sharp interface tracking
        // Each region corresponds to a slab's gas properties
        this.regions = slabs.map(slab => ({
            gamma: slab.gamma,
            mw: slab.mw,
            gasId: slab.gasId
        }));
        
        // Initialize Lagrangian tracers at slab interfaces
        currentPos = 0;
        for (let i = 0; i < slabs.length; i++) {
            currentPos += slabs[i].length;
            if (i < slabs.length - 1) {  // Don't add tracer at end
                this.tracers.push({
                    x: currentPos,
                    trajectory: [{t: 0, x: currentPos}]
                });
            }
        }
        
        // Update primitive variables
        this.updatePrimitives();
        
        // Store initial x-t data
        this.storeXTData();
    }
    
    /**
     * Convert conservative to primitive variables
     */
    updatePrimitives() {
        for (let i = 0; i < this.nx; i++) {
            const idx = i * 3;
            const rho = this.U[idx];
            const rhou = this.U[idx + 1];
            const E = this.U[idx + 2];
            const gamma = this.gamma[i];
            
            this.rho[i] = rho;
            this.u[i] = rhou / rho;
            
            // p = (gamma - 1) * (E - 0.5 * rho * u^2)
            const u = this.u[i];
            this.p[i] = (gamma - 1) * (E - 0.5 * rho * u * u);
            
            // T = p / (rho * R)
            const R = this.Ru / this.mw[i];
            this.T[i] = this.p[i] / (rho * R);
        }
    }
    
    /**
     * Compute HLLC flux at cell interface
     */
    hllcFlux(UL, UR, gammaL, gammaR, fluxOut, idx) {
        // Left state
        const rhoL = UL[0];
        const uL = UL[1] / rhoL;
        const EL = UL[2];
        const pL = (gammaL - 1) * (EL - 0.5 * rhoL * uL * uL);
        const aL = Math.sqrt(gammaL * pL / rhoL);
        const HL = (EL + pL) / rhoL;
        
        // Right state
        const rhoR = UR[0];
        const uR = UR[1] / rhoR;
        const ER = UR[2];
        const pR = (gammaR - 1) * (ER - 0.5 * rhoR * uR * uR);
        const aR = Math.sqrt(gammaR * pR / rhoR);
        const HR = (ER + pR) / rhoR;
        
        // Wave speed estimates (Davis)
        const SL = Math.min(uL - aL, uR - aR);
        const SR = Math.max(uL + aL, uR + aR);
        
        // Middle wave speed
        const SStar = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) /
                      (rhoL * (SL - uL) - rhoR * (SR - uR));
        
        // Compute flux based on wave structure
        if (SL >= 0) {
            // Left flux
            fluxOut[idx] = rhoL * uL;
            fluxOut[idx + 1] = rhoL * uL * uL + pL;
            fluxOut[idx + 2] = uL * (EL + pL);
        } else if (SStar >= 0) {
            // Left star flux
            const pStar = pL + rhoL * (SL - uL) * (SStar - uL);
            const rhoStarL = rhoL * (SL - uL) / (SL - SStar);
            const EStarL = rhoStarL * (EL / rhoL + (SStar - uL) * (SStar + pL / (rhoL * (SL - uL))));
            
            fluxOut[idx] = rhoStarL * SStar;
            fluxOut[idx + 1] = rhoStarL * SStar * SStar + pStar;
            fluxOut[idx + 2] = SStar * (EStarL + pStar);
        } else if (SR >= 0) {
            // Right star flux
            const pStar = pR + rhoR * (SR - uR) * (SStar - uR);
            const rhoStarR = rhoR * (SR - uR) / (SR - SStar);
            const EStarR = rhoStarR * (ER / rhoR + (SStar - uR) * (SStar + pR / (rhoR * (SR - uR))));
            
            fluxOut[idx] = rhoStarR * SStar;
            fluxOut[idx + 1] = rhoStarR * SStar * SStar + pStar;
            fluxOut[idx + 2] = SStar * (EStarR + pStar);
        } else {
            // Right flux
            fluxOut[idx] = rhoR * uR;
            fluxOut[idx + 1] = rhoR * uR * uR + pR;
            fluxOut[idx + 2] = uR * (ER + pR);
        }
        
        return Math.max(Math.abs(SL), Math.abs(SR));
    }
    
    /**
     * Compute fluxes at all cell interfaces
     */
    computeFluxes() {
        let maxWaveSpeed = 0;
        
        for (let i = 0; i <= this.nx; i++) {
            let UL = new Float64Array(3);
            let UR = new Float64Array(3);
            let gammaL, gammaR;
            
            if (i === 0) {
                // Left boundary - reflective
                const idx = 0;
                UL[0] = this.U[idx];
                UL[1] = -this.U[idx + 1];  // Reflect velocity
                UL[2] = this.U[idx + 2];
                UR[0] = this.U[idx];
                UR[1] = this.U[idx + 1];
                UR[2] = this.U[idx + 2];
                gammaL = this.gamma[0];
                gammaR = this.gamma[0];
            } else if (i === this.nx) {
                // Right boundary - reflective
                const idx = (this.nx - 1) * 3;
                UL[0] = this.U[idx];
                UL[1] = this.U[idx + 1];
                UL[2] = this.U[idx + 2];
                UR[0] = this.U[idx];
                UR[1] = -this.U[idx + 1];  // Reflect velocity
                UR[2] = this.U[idx + 2];
                gammaL = this.gamma[this.nx - 1];
                gammaR = this.gamma[this.nx - 1];
            } else {
                // Interior interface
                const idxL = (i - 1) * 3;
                const idxR = i * 3;
                UL[0] = this.U[idxL];
                UL[1] = this.U[idxL + 1];
                UL[2] = this.U[idxL + 2];
                UR[0] = this.U[idxR];
                UR[1] = this.U[idxR + 1];
                UR[2] = this.U[idxR + 2];
                gammaL = this.gamma[i - 1];
                gammaR = this.gamma[i];
            }
            
            const waveSpeed = this.hllcFlux(UL, UR, gammaL, gammaR, this.flux, i * 3);
            maxWaveSpeed = Math.max(maxWaveSpeed, waveSpeed);
        }
        
        return maxWaveSpeed;
    }
    
    /**
     * Update conservative variables using computed fluxes
     */
    updateConservative(U, Uout, dt) {
        for (let i = 0; i < this.nx; i++) {
            const idx = i * 3;
            const fluxL = (i) * 3;
            const fluxR = (i + 1) * 3;
            
            for (let k = 0; k < 3; k++) {
                Uout[idx + k] = U[idx + k] - dt / this.dx * (this.flux[fluxR + k] - this.flux[fluxL + k]);
            }
        }
    }
    
    /**
     * Take one time step using RK2
     */
    step() {
        // Stage 1: compute flux and update
        const maxWaveSpeed = this.computeFluxes();
        
        // Adaptive time step based on CFL condition
        this.dt = this.cfl * this.dx / maxWaveSpeed;
        
        // Don't overshoot final time
        if (this.t + this.dt > this.finalTime) {
            this.dt = this.finalTime - this.t;
        }
        
        // RK Stage 1: U* = U^n + dt * L(U^n)
        this.updateConservative(this.U, this.Uk, this.dt);
        
        // Stage 2: compute flux from U*
        // Copy Uk to U temporarily for flux computation
        for (let i = 0; i < this.nx * 3; i++) {
            this.Un[i] = this.U[i];
            this.U[i] = this.Uk[i];
        }
        
        this.computeFluxes();
        
        // RK Stage 2: U^{n+1} = 0.5 * U^n + 0.5 * (U* + dt * L(U*))
        this.updateConservative(this.Uk, this.Uk, this.dt);
        
        for (let i = 0; i < this.nx * 3; i++) {
            this.U[i] = 0.5 * this.Un[i] + 0.5 * this.Uk[i];
        }
        
        // Update time and step counter
        this.t += this.dt;
        this.timeSteps++;
        
        // Update primitive variables
        this.updatePrimitives();
        
        // Update Lagrangian tracers
        this.updateTracers();
        
        // Update gas properties based on interface positions (sharp interface tracking)
        this.updateGasProperties();
        
        // Store x-t data if needed
        if (this.t >= this.nextXTStore) {
            this.storeXTData();
            this.nextXTStore += this.xtInterval;
        }
    }
    
    /**
     * Update gas properties in each cell based on interface positions
     * Uses sharp interface tracking - cells inherit properties from their region
     */
    updateGasProperties() {
        // Sort tracer positions for efficient region determination
        const tracerPositions = this.tracers.map(t => t.x).sort((a, b) => a - b);
        
        for (let i = 0; i < this.nx; i++) {
            const cellCenter = this.x[i];
            
            // Determine which region this cell is in
            let regionIndex = 0;
            for (let j = 0; j < tracerPositions.length; j++) {
                if (cellCenter > tracerPositions[j]) {
                    regionIndex = j + 1;
                } else {
                    break;
                }
            }
            
            // Clamp region index to valid range
            regionIndex = Math.max(0, Math.min(this.regions.length - 1, regionIndex));
            
            // Update cell properties from the appropriate region
            const region = this.regions[regionIndex];
            this.gamma[i] = region.gamma;
            this.mw[i] = region.mw;
            this.gasId[i] = region.gasId;
        }
    }
    
    /**
     * Update positions of Lagrangian tracers
     */
    updateTracers() {
        for (let tracer of this.tracers) {
            // Find cell containing tracer
            const cellIdx = Math.floor(tracer.x / this.dx);
            const clampedIdx = Math.max(0, Math.min(this.nx - 1, cellIdx));
            
            // Get velocity at tracer location
            const velocity = this.u[clampedIdx];
            
            // Update position using forward Euler (good enough for tracers)
            tracer.x += velocity * this.dt;
            
            // Clamp to domain
            tracer.x = Math.max(0, Math.min(this.x[this.nx - 1] + 0.5 * this.dx, tracer.x));
            
            // Store trajectory point
            tracer.trajectory.push({t: this.t, x: tracer.x});
        }
    }
    
    /**
     * Store current pressure field for x-t diagram
     */
    storeXTData() {
        const snapshot = {
            t: this.t,
            p: new Float64Array(this.p)  // Copy pressure array
        };
        this.xtData.push(snapshot);
    }
    
    /**
     * Run simulation to completion
     */
    run(progressCallback) {
        const startTime = Date.now();
        let lastProgressUpdate = startTime;
        
        while (this.t < this.finalTime) {
            this.step();
            
            // Update progress every 100ms
            const now = Date.now();
            if (progressCallback && now - lastProgressUpdate > 100) {
                progressCallback(this.t / this.finalTime);
                lastProgressUpdate = now;
            }
        }
        
        // Final progress update
        if (progressCallback) {
            progressCallback(1.0);
        }
        
        const endTime = Date.now();
        console.log(`Simulation completed in ${(endTime - startTime) / 1000} seconds`);
        console.log(`Time steps: ${this.timeSteps}`);
        console.log(`Final time: ${this.t} s`);
        console.log(`X-T snapshots: ${this.xtData.length}`);
    }
    
    /**
     * Get results for export
     */
    getResults() {
        return {
            x: Array.from(this.x),
            t: this.t,
            rho: Array.from(this.rho),
            u: Array.from(this.u),
            p: Array.from(this.p),
            T: Array.from(this.T),
            xtData: this.xtData,
            tracers: this.tracers,
            timeSteps: this.timeSteps
        };
    }
}
