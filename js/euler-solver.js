/**
 * 1D Compressible Euler Solver with HLLC Flux
 * For multicomponent shock tube simulations
 */

class EulerSolver {
    constructor(config) {
        this.nx = config.nx || 500;  // Number of grid cells
        this.cfl = config.cfl || 0.4;  // CFL number for stability
        this.finalTime = config.finalTime || 0.02;  // Final simulation time (seconds)
        
        // Interface tracking method: 'sharp', 'ghost', or 'mixed'
        this.interfaceMethod = config.interfaceMethod || 'sharp';
        
        // Create time integrator instance
        this.integrator = IntegratorFactory.create(config.integrator || 'RK2');

        // Initialize arrays using typed arrays for performance
        this.x = new Float64Array(this.nx);  // Cell centers
        this.dx = 0;  // Grid spacing (set during initialization)
        
        // Initialize conservative variable buffers based on integrator requirements
        // All integrators now use the same 3-buffer structure: U, Un, Uk
        this.U = new Float64Array(this.nx * 3);   // Current state
        this.Un = new Float64Array(this.nx * 3);  // Storage buffer
        this.Uk = new Float64Array(this.nx * 3);  // Work buffer
        
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
        
        // Material fraction tracking for mixed-cell method
        this.numMaterials = 0;  // Set during initialization
        this.materialFractions = null;  // Allocated if interfaceMethod === 'mixed'
        this.materialProperties = [];  // Array of {gamma, mw, gasId} for each material
        this.materialFlux = null;  // Flux for material fractions
        
        // Ghost fluid method tracking
        this.interfaceCells = new Set();  // Cells adjacent to material interfaces
        
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
        
        // Initialize mixed-cell method if selected
        if (this.interfaceMethod === 'mixed') {
            this.numMaterials = slabs.length;
            this.materialFractions = new Float64Array(this.nx * this.numMaterials);
            this.materialFlux = new Float64Array((this.nx + 1) * this.numMaterials);
            
            // Store material properties
            this.materialProperties = slabs.map(slab => ({
                gamma: slab.gamma,
                mw: slab.mw,
                gasId: slab.gasId
            }));
            
            // Initialize material fractions based on slab positions
            currentPos = 0;
            slabIndex = 0;
            for (let i = 0; i < this.nx; i++) {
                const cellLeft = i * this.dx;
                const cellRight = (i + 1) * this.dx;
                
                // Find which slab(s) this cell overlaps
                let tempPos = 0;
                for (let s = 0; s < slabs.length; s++) {
                    const slabLeft = tempPos;
                    const slabRight = tempPos + slabs[s].length;
                    
                    // Calculate overlap
                    const overlapLeft = Math.max(cellLeft, slabLeft);
                    const overlapRight = Math.min(cellRight, slabRight);
                    const overlap = Math.max(0, overlapRight - overlapLeft);
                    
                    this.materialFractions[i * this.numMaterials + s] = overlap / this.dx;
                    
                    tempPos += slabs[s].length;
                }
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
     * Compute mixture properties from material fractions
     * Uses thermodynamically consistent mixing rules
     */
    computeMixtureProperties(cellIdx) {
        let gamma_mix_inv = 0;
        let mw_mix = 0;
        let totalFraction = 0;
        let dominantMaterial = 0;
        let maxFraction = 0;
        
        for (let m = 0; m < this.numMaterials; m++) {
            const alpha = this.materialFractions[cellIdx * this.numMaterials + m];
            const gamma = this.materialProperties[m].gamma;
            const mw = this.materialProperties[m].mw;
            
            // Pressure-weighted average for gamma (thermodynamically consistent)
            gamma_mix_inv += alpha / (gamma - 1);
            
            // Mass-weighted average for molecular weight
            mw_mix += alpha * mw;
            
            totalFraction += alpha;
            
            // Track dominant material for gas ID
            if (alpha > maxFraction) {
                maxFraction = alpha;
                dominantMaterial = m;
            }
        }
        
        // Normalize (should be ~1 but numerical errors possible)
        if (totalFraction > 1e-10) {
            gamma_mix_inv /= totalFraction;
            mw_mix /= totalFraction;
        } else {
            // Fallback to first material if fractions sum to zero
            gamma_mix_inv = 1.0 / (this.materialProperties[0].gamma - 1);
            mw_mix = this.materialProperties[0].mw;
            dominantMaterial = 0;
        }
        
        return {
            gamma: 1 + 1 / gamma_mix_inv,
            mw: mw_mix,
            gasId: this.materialProperties[dominantMaterial].gasId
        };
    }
    
    /**
     * Identify cells adjacent to material interfaces for ghost fluid method
     */
    identifyInterfaceCells() {
        this.interfaceCells.clear();
        
        // Sort tracer positions
        const tracerPositions = this.tracers.map(t => t.x).sort((a, b) => a - b);
        
        for (let i = 0; i < this.nx; i++) {
            const cellLeft = i * this.dx;
            const cellRight = (i + 1) * this.dx;
            
            // Check if any tracer is within or near this cell
            for (const tracerX of tracerPositions) {
                const dist = Math.min(
                    Math.abs(tracerX - cellLeft),
                    Math.abs(tracerX - cellRight),
                    Math.abs(tracerX - this.x[i])
                );
                
                // Mark as interface cell if tracer is within 1.5 cell widths
                if (dist < 1.5 * this.dx) {
                    this.interfaceCells.add(i);
                    break;
                }
            }
        }
    }
    
    /**
     * Update gas properties in each cell based on interface tracking method
     */
    updateGasProperties() {
        if (this.interfaceMethod === 'mixed') {
            // Mixed-cell method: compute mixture properties from volume fractions
            for (let i = 0; i < this.nx; i++) {
                const mixProps = this.computeMixtureProperties(i);
                this.gamma[i] = mixProps.gamma;
                this.mw[i] = mixProps.mw;
                this.gasId[i] = mixProps.gasId;
            }
        } else if (this.interfaceMethod === 'ghost') {
            // Ghost fluid method: sharp interfaces with special treatment near boundaries
            // First, identify interface cells
            this.identifyInterfaceCells();
            
            // Sort tracer positions for region determination
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
            
            // For ghost fluid method, we could apply special treatment to interface cells here
            // For now, the sharp assignment is sufficient but can be enhanced later
        } else {
            // Sharp method: original behavior - cells inherit properties from their region
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
    }
    
    /**
     * Advect material fractions using upwind scheme
     * Solves ∂α/∂t + u·∂α/∂x = 0 for each material
     */
    advectMaterialFractions(dt) {
        if (this.interfaceMethod !== 'mixed' || !this.materialFractions) {
            return;
        }
        
        // Compute fluxes for material fractions at cell interfaces
        for (let i = 0; i <= this.nx; i++) {
            for (let m = 0; m < this.numMaterials; m++) {
                const fluxIdx = i * this.numMaterials + m;
                
                if (i === 0 || i === this.nx) {
                    // Boundary: zero flux (reflective)
                    this.materialFlux[fluxIdx] = 0;
                } else {
                    // Interior interface: upwind scheme
                    // Interface velocity is average of adjacent cells
                    const uInterface = 0.5 * (this.u[i - 1] + this.u[i]);
                    
                    if (uInterface > 0) {
                        // Upwind from left
                        const alphaL = this.materialFractions[(i - 1) * this.numMaterials + m];
                        this.materialFlux[fluxIdx] = uInterface * alphaL;
                    } else {
                        // Upwind from right
                        const alphaR = this.materialFractions[i * this.numMaterials + m];
                        this.materialFlux[fluxIdx] = uInterface * alphaR;
                    }
                }
            }
        }
        
        // Update material fractions using computed fluxes
        const alphaNew = new Float64Array(this.nx * this.numMaterials);
        
        for (let i = 0; i < this.nx; i++) {
            for (let m = 0; m < this.numMaterials; m++) {
                const idx = i * this.numMaterials + m;
                const fluxL = i * this.numMaterials + m;
                const fluxR = (i + 1) * this.numMaterials + m;
                
                alphaNew[idx] = this.materialFractions[idx] - 
                               dt / this.dx * (this.materialFlux[fluxR] - this.materialFlux[fluxL]);
            }
            
            // Normalize fractions to ensure they sum to 1
            let totalAlpha = 0;
            for (let m = 0; m < this.numMaterials; m++) {
                const idx = i * this.numMaterials + m;
                alphaNew[idx] = Math.max(0, Math.min(1, alphaNew[idx]));  // Clamp to [0,1]
                totalAlpha += alphaNew[idx];
            }
            
            // Renormalize
            if (totalAlpha > 1e-10) {
                for (let m = 0; m < this.numMaterials; m++) {
                    const idx = i * this.numMaterials + m;
                    alphaNew[idx] /= totalAlpha;
                }
            } else {
                // If all fractions are zero, default to equal distribution
                for (let m = 0; m < this.numMaterials; m++) {
                    const idx = i * this.numMaterials + m;
                    alphaNew[idx] = 1.0 / this.numMaterials;
                }
            }
        }
        
        // Copy updated fractions back
        this.materialFractions.set(alphaNew);
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
     * Store current state for x-t diagram
     * Includes all primitive variables and material properties
     */
    storeXTData() {
        const snapshot = {
            t: this.t,
            rho: new Float64Array(this.rho),    // Density
            u: new Float64Array(this.u),        // Velocity
            p: new Float64Array(this.p),        // Pressure
            T: new Float64Array(this.T),        // Temperature
            gamma: new Float64Array(this.gamma), // Specific heat ratio
            mw: new Float64Array(this.mw),      // Molecular weight
            gasId: Array.from(this.gasId)       // Gas identifier
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
            // Delegate time stepping to the integrator
            this.integrator.step(this);
            
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
