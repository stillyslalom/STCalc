/**
 * 1D Compressible Euler Solver with HLLC Flux
 * For multicomponent shock tube simulations with ideal or stiffened gases
 */

class EulerSolver {
    constructor(config) {
        this.nx = config.nx || 500;  // Number of grid cells
        this.cfl = config.cfl || 0.4;  // CFL number for stability
        this.finalTime = config.finalTime || 0.02;  // Final simulation time (seconds)

        // Interface tracking method: 'sharp' or 'mixed'
        this.interfaceMethod = config.interfaceMethod || 'sharp';

        // Create time integrator instance
        this.integrator = IntegratorFactory.create(config.integrator || 'RK2');

        this.x = new Float64Array(this.nx);  // Cell centers
        this.dx = 0;  // Grid spacing (set during initialization)

        // State buffers are allocated in initialize(), once the number of
        // materials (and therefore variables per cell) is known.
        // Layout per cell: [rho, rho*u, E, alpha_0 ... alpha_{M-1}]
        // (volume fractions are only carried by the mixed-cell method), plus
        // partial densities alpha_m * rho_m when any material is tabulated
        this.nvar = 3;
        this.U = null;   // Current state
        this.Un = null;  // Storage buffer
        this.Uk = null;  // Work buffer

        // Primitive variables for output
        this.rho = new Float64Array(this.nx);
        this.u = new Float64Array(this.nx);
        this.p = new Float64Array(this.nx);
        this.T = new Float64Array(this.nx);

        // Effective equation-of-state properties per cell
        this.gamma = new Float64Array(this.nx);
        this.mw = new Float64Array(this.nx);    // Molecular weight (volume-fraction weighted in mixed cells)
        this.pinf = new Float64Array(this.nx);  // Stiffening pressure
        this.gasId = new Array(this.nx);        // Gas identifier for each cell

        // Stiffened-gas closure used by the flux: E = G*p + P + 0.5*rho*u^2
        // with G = 1/(gamma-1) and P = gamma*pinf/(gamma-1)
        this.eosG = new Float64Array(this.nx);
        this.eosP = new Float64Array(this.nx);

        // Cell pressure and rho*a^2 for the current flux evaluation
        this.pc = new Float64Array(this.nx);
        this.rc2 = new Float64Array(this.nx);

        // True when any material has a tabulated (temperature-dependent)
        // equation of state; partial densities are then carried per material
        this.tabulated = false;

        // Face arrays
        this.flux = null;
        this.faceVelocity = new Float64Array(this.nx + 1);

        // Time tracking
        this.t = 0;
        this.dt = 0;
        this.timeSteps = 0;

        // X-T diagram storage
        this.xtData = [];
        this.xtInterval = config.xtInterval || 0.0001;  // Store every 0.1 ms
        this.nextXTStore = 0;

        // Lagrangian tracers
        this.tracers = [];

        // Regions for sharp interface tracking
        this.regions = [];

        // Materials for the mixed-cell method
        this.numMaterials = 0;
        this.materialFractions = null;  // Output copy of the volume fractions
        this.materialProperties = [];   // {gamma, mw, pinf, gasId} per material

        // Constants
        this.Ru = 8314.51;  // Universal gas constant J/(kmol·K)
    }

    /**
     * Initialize the shock tube with gas slabs
     * @param {Array} slabs - Array of {gasId, gamma, mw, pinf, pressure, temperature, length, velocity}.
     *   A slab may also carry thermo: a composition for makeThermo() (thermo.js);
     *   its specific heats then depend on temperature.
     */
    initialize(slabs) {
        const mixed = this.interfaceMethod === 'mixed';
        const materials = slabs.map(slab => ({
            gamma: slab.gamma,
            mw: slab.mw,
            pinf: slab.pinf || 0,
            gasId: slab.gasId,
            G: 1 / (slab.gamma - 1),
            P: slab.gamma * (slab.pinf || 0) / (slab.gamma - 1),
            R: this.Ru / slab.mw,
            table: null
        }));

        if (!mixed && materials.some(m => m.pinf !== materials[0].pinf)) {
            throw new Error('Slabs with different stiffening pressures (e.g. gas next to liquid) require the Mixed-Cell interface method.');
        }

        slabs.forEach((slab, s) => {
            if (!slab.thermo || materials[s].pinf > 0) return;
            if (!mixed) throw new Error('Temperature-dependent gas properties require the Mixed-Cell interface method.');
            const th = makeThermo(slab.thermo);
            materials[s].R = th.R;
            materials[s].mw = th.mw;
            materials[s].table = EulerSolver.energyTable(th);
            this.tabulated = true;
        });

        this.materialProperties = materials;
        this.regions = materials;
        this.numMaterials = mixed ? slabs.length : 0;
        this.nvar = 3 + this.numMaterials * (this.tabulated ? 2 : 1);

        const n = this.nx * this.nvar;
        this.U = new Float64Array(n);
        this.Un = new Float64Array(n);
        this.Uk = new Float64Array(n);
        this.flux = new Float64Array((this.nx + 1) * this.nvar);
        if (mixed) {
            this.materialFractions = new Float64Array(this.nx * this.numMaterials);
        }

        const totalLength = slabs.reduce((sum, slab) => sum + slab.length, 0);
        this.dx = totalLength / this.nx;
        for (let i = 0; i < this.nx; i++) {
            this.x[i] = (i + 0.5) * this.dx;
        }

        // Slab edges
        const edges = [0];
        for (const slab of slabs) edges.push(edges[edges.length - 1] + slab.length);

        for (let i = 0; i < this.nx; i++) {
            const idx = i * this.nvar;
            const cellLeft = i * this.dx;
            const cellRight = (i + 1) * this.dx;

            if (mixed) {
                // Volume-average every slab overlapping this cell so that
                // density, energy and volume fractions are mutually consistent
                let rho = 0, mom = 0, internal = 0, total = 0;
                for (let s = 0; s < slabs.length; s++) {
                    const overlap = Math.max(0, Math.min(cellRight, edges[s + 1]) - Math.max(cellLeft, edges[s]));
                    const alpha = overlap / this.dx;
                    this.U[idx + 3 + s] = alpha;
                    if (alpha <= 0) continue;
                    const mat = materials[s];
                    const rhoS = (slabs[s].pressure + mat.pinf) / (mat.R * slabs[s].temperature);
                    const uS = slabs[s].velocity || 0;
                    rho += alpha * rhoS;
                    mom += alpha * rhoS * uS;
                    internal += alpha * (mat.table
                        ? rhoS * EulerSolver.tableEnergy(mat.table, slabs[s].temperature)
                        : mat.G * slabs[s].pressure + mat.P);
                    total += alpha;
                    if (this.tabulated) this.U[idx + 3 + slabs.length + s] = alpha * rhoS;
                }
                // Guard against round-off at the domain end
                for (let s = 0; s < slabs.length; s++) {
                    this.U[idx + 3 + s] /= total;
                    if (this.tabulated) this.U[idx + 3 + slabs.length + s] /= total;
                }
                this.pc[i] = slabs.reduce((sum, slab, s) => sum + this.U[idx + 3 + s] * slab.pressure, 0);
                this.U[idx] = rho / total;
                this.U[idx + 1] = mom / total;
                this.U[idx + 2] = internal / total + 0.5 * mom * mom / (rho * total);
            } else {
                let s = 0;
                while (s < slabs.length - 1 && this.x[i] >= edges[s + 1]) s++;
                const mat = materials[s];
                const rho = (slabs[s].pressure + mat.pinf) / (mat.R * slabs[s].temperature);
                const u = slabs[s].velocity || 0;
                this.U[idx] = rho;
                this.U[idx + 1] = rho * u;
                this.U[idx + 2] = mat.G * slabs[s].pressure + mat.P + 0.5 * rho * u * u;
            }
        }

        // Lagrangian tracers at slab interfaces
        for (let s = 1; s < slabs.length; s++) {
            this.tracers.push({
                x: edges[s],
                trajectory: [{t: 0, x: edges[s], v: 0}]
            });
        }

        this.updateGasProperties();
        this.updatePrimitives();
        this.storeXTData();
    }

    /**
     * Refresh the per-cell equation of state (eosG, eosP) from a state vector.
     * Mixed-cell: volume-fraction average of 1/(gamma-1) and gamma*pinf/(gamma-1),
     * which keeps pressure and velocity uniform across material interfaces.
     * Sharp: properties of the region between Lagrangian tracers.
     */
    updateEOS(U) {
        const nvar = this.nvar;
        if (this.tabulated) {
            this.closeTabulated(U);
        } else if (this.interfaceMethod === 'mixed') {
            const M = this.numMaterials;
            const mats = this.materialProperties;
            for (let i = 0; i < this.nx; i++) {
                const base = i * nvar + 3;
                let G = 0, P = 0, total = 0;
                for (let m = 0; m < M; m++) {
                    const alpha = Math.min(1, Math.max(0, U[base + m]));
                    G += alpha * mats[m].G;
                    P += alpha * mats[m].P;
                    total += alpha;
                }
                if (total > 1e-12) {
                    this.eosG[i] = G / total;
                    this.eosP[i] = P / total;
                } else {
                    this.eosG[i] = mats[0].G;
                    this.eosP[i] = mats[0].P;
                }
            }
        } else {
            const tracerPositions = this.tracers.map(t => t.x).sort((a, b) => a - b);
            let regionIndex = 0;
            for (let i = 0; i < this.nx; i++) {
                while (regionIndex < tracerPositions.length && this.x[i] > tracerPositions[regionIndex]) {
                    regionIndex++;
                }
                const region = this.regions[Math.min(regionIndex, this.regions.length - 1)];
                this.eosG[i] = region.G;
                this.eosP[i] = region.P;
            }
        }
    }

    /**
     * Specific internal energy e(T) of a thermally perfect gas, tabulated on a
     * uniform temperature grid and interpolated linearly. e(0) = 0, as for a
     * constant-gamma gas. A piecewise-linear e makes the pressure solve cheap
     * and exactly invertible.
     */
    static energyTable(th) {
        const dT = 2, n = 12001;                  // 0 to 24000 K
        const e = new Float64Array(n);
        // Integrate cv with the midpoint rule; cp is held constant below the fitted range
        for (let j = 1; j < n; j++) {
            e[j] = e[j - 1] + (th.cp((j - 0.5) * dT) - th.R) * dT;
        }
        return { dT, n, e };
    }

    static tableEnergy(table, T) {
        const x = T / table.dT;
        const j = Math.min(table.n - 2, Math.max(0, Math.floor(x)));
        return table.e[j] + (table.e[j + 1] - table.e[j]) * (x - j);
    }

    /** cv(T), the slope of the tabulated energy. */
    static tableCv(table, T) {
        const j = Math.min(table.n - 2, Math.max(0, Math.floor(T / table.dT)));
        return (table.e[j + 1] - table.e[j]) / table.dT;
    }

    /**
     * Pressure closure when materials have tabulated equations of state.
     * Each material k occupies a volume fraction alpha_k at its own density
     * rho_k and temperature, and all share one pressure p, found from
     *
     *     sum_k alpha_k * eps_k(p, rho_k) = rho * e
     *
     * where eps_k is internal energy per unit volume: G_k p + P_k for a
     * constant-gamma or stiffened material, rho_k e_k(p / (rho_k R_k)) for a
     * tabulated one. Fills pc, rc2 and the effective eosG, eosP.
     */
    closeTabulated(U) {
        const nvar = this.nvar;
        const M = this.numMaterials;
        const mats = this.materialProperties;
        const tiny = 1e-10;

        for (let i = 0; i < this.nx; i++) {
            const idx = i * nvar;
            const rho = Math.max(U[idx], tiny);
            const u = U[idx + 1] / rho;
            const target = U[idx + 2] - 0.5 * rho * u * u;

            // Lowest admissible pressure for the materials present
            let floor = -Infinity;
            for (let m = 0; m < M; m++) {
                if (U[idx + 3 + m] > 1e-12) floor = Math.max(floor, -mats[m].pinf);
            }
            floor += tiny;

            let p = Math.max(this.pc[i], floor);
            let xi = 1, stiff = 0;
            for (let iter = 0; iter < 40; iter++) {
                let energy = 0;
                xi = 0;
                stiff = 0;
                for (let m = 0; m < M; m++) {
                    const alpha = U[idx + 3 + m];
                    if (!(alpha > 1e-12)) continue;
                    const mat = mats[m];
                    const mass = U[idx + 3 + M + m];
                    if (mat.table && mass > 0) {
                        const rhoK = mass / alpha;
                        const T = p / (rhoK * mat.R);
                        energy += mass * EulerSolver.tableEnergy(mat.table, T);
                        xi += alpha * EulerSolver.tableCv(mat.table, T) / mat.R;
                    } else {
                        energy += alpha * (mat.G * p + mat.P);
                        xi += alpha * mat.G;
                        stiff += alpha * mat.P;
                    }
                }
                if (!(xi > 0)) break;
                let next = p - (energy - target) / xi;
                if (!(next > floor)) next = 0.5 * (p + floor);
                const done = Math.abs(next - p) <= 1e-12 * Math.abs(p) + 1e-9;
                p = next;
                if (done) break;
            }

            this.pc[i] = p;
            this.eosG[i] = xi;
            this.eosP[i] = stiff;
            // Mixture sound speed: rho a^2 = [p * sum alpha (xi_k + 1) + sum alpha P_k] / xi
            let alphaSum = 0;
            for (let m = 0; m < M; m++) if (U[idx + 3 + m] > 1e-12) alphaSum += U[idx + 3 + m];
            this.rc2[i] = Math.max((p * (xi + alphaSum) + stiff) / xi, tiny);
        }
    }

    /**
     * Convert conservative to primitive variables
     */
    updatePrimitives() {
        const minDensity = 1e-10;  // Minimum density threshold for numerical stability
        const nvar = this.nvar;
        const mixed = this.interfaceMethod === 'mixed';

        for (let i = 0; i < this.nx; i++) {
            const idx = i * nvar;
            const rho = Math.max(this.U[idx], minDensity);
            const u = this.U[idx + 1] / rho;
            const G = this.eosG[i];
            const P = this.eosP[i];
            const pinf = P / (G + 1);

            // p = (E - 0.5*rho*u^2 - P) / G, floored so that p + pinf > 0
            let p = (this.U[idx + 2] - 0.5 * rho * u * u - P) / G;
            const pFloor = -pinf + minDensity;
            if (!(p > pFloor)) p = pFloor;
            if (this.tabulated) p = this.pc[i];   // from closeTabulated

            this.rho[i] = rho;
            this.u[i] = u;
            this.p[i] = p;
            this.gamma[i] = 1 + 1 / G;
            this.pinf[i] = pinf;

            if (mixed) {
                // Mixture temperature assuming the materials in a cell share p and T
                const M = this.numMaterials;
                const mats = this.materialProperties;
                let sum = 0, mw = 0, total = 0, dominant = 0, maxAlpha = -1;
                for (let m = 0; m < M; m++) {
                    const alpha = Math.min(1, Math.max(0, this.U[idx + 3 + m]));
                    this.materialFractions[i * M + m] = alpha;
                    // alpha_k * rho_k * T_k: each material is at its own temperature
                    sum += alpha * (p + mats[m].pinf) / mats[m].R;
                    mw += alpha * mats[m].mw;
                    total += alpha;
                    if (alpha > maxAlpha) {
                        maxAlpha = alpha;
                        dominant = m;
                    }
                }
                this.T[i] = sum / (total * rho);
                this.mw[i] = mw / total;
                this.gasId[i] = mats[dominant].gasId;
            } else {
                // Region lookup mirrors updateEOS
                const region = this.regionOfCell(i);
                this.T[i] = (p + pinf) / (rho * region.R);
                this.mw[i] = region.mw;
                this.gasId[i] = region.gasId;
            }
        }
    }

    /** Region (slab material) containing cell i for the sharp method. */
    regionOfCell(i) {
        let regionIndex = 0;
        for (let j = 0; j < this.tracers.length; j++) {
            if (this.x[i] > this.tracers[j].x) regionIndex = j + 1;
        }
        return this.regions[Math.min(regionIndex, this.regions.length - 1)];
    }

    /**
     * HLLC flux between a left and a right state, each given by its
     * conservative variables, pressure and rho*a^2.
     * Writes the conservative flux to this.flux at offset fluxIdx and returns
     * the largest wave speed. The velocity of the Riemann solution at the face
     * is left in this._uFace and the upwind side (0 = left, 1 = right) of the
     * contact in this._upwind; both are needed to advect volume fractions.
     */
    hllcFlux(rhoL, momL, EL, pL, rc2L, rhoR, momR, ER, pR, rc2R, fluxIdx) {
        const tiny = 1e-10;
        const flux = this.flux;

        rhoL = Math.max(rhoL, tiny);
        rhoR = Math.max(rhoR, tiny);
        const uL = momL / rhoL;
        const uR = momR / rhoR;
        const aL = Math.sqrt(rc2L / rhoL);
        const aR = Math.sqrt(rc2R / rhoR);

        // Wave speed estimates (Davis)
        const SL = Math.min(uL - aL, uR - aR);
        const SR = Math.max(uL + aL, uR + aR);

        // Contact wave speed
        const SStar = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) /
                      (rhoL * (SL - uL) - rhoR * (SR - uR));

        if (SL >= 0) {
            flux[fluxIdx] = rhoL * uL;
            flux[fluxIdx + 1] = rhoL * uL * uL + pL;
            flux[fluxIdx + 2] = uL * (EL + pL);
            this._uFace = uL;
        } else if (SStar >= 0) {
            const pStar = pL + rhoL * (SL - uL) * (SStar - uL);
            const rhoStarL = rhoL * (SL - uL) / (SL - SStar);
            const EStarL = rhoStarL * (EL / rhoL + (SStar - uL) * (SStar + pL / (rhoL * (SL - uL))));

            flux[fluxIdx] = rhoStarL * SStar;
            flux[fluxIdx + 1] = rhoStarL * SStar * SStar + pStar;
            flux[fluxIdx + 2] = SStar * (EStarL + pStar);
            this._uFace = SStar;
        } else if (SR >= 0) {
            const pStar = pR + rhoR * (SR - uR) * (SStar - uR);
            const rhoStarR = rhoR * (SR - uR) / (SR - SStar);
            const EStarR = rhoStarR * (ER / rhoR + (SStar - uR) * (SStar + pR / (rhoR * (SR - uR))));

            flux[fluxIdx] = rhoStarR * SStar;
            flux[fluxIdx + 1] = rhoStarR * SStar * SStar + pStar;
            flux[fluxIdx + 2] = SStar * (EStarR + pStar);
            this._uFace = SStar;
        } else {
            flux[fluxIdx] = rhoR * uR;
            flux[fluxIdx + 1] = rhoR * uR * uR + pR;
            flux[fluxIdx + 2] = uR * (ER + pR);
            this._uFace = uR;
        }
        this._upwind = SStar >= 0 ? 0 : 1;

        return Math.max(Math.abs(SL), Math.abs(SR));
    }

    /**
     * Compute fluxes at all cell interfaces from the current state this.U.
     * Volume-fraction entries of the flux hold alpha_upwind * u_face; the
     * non-conservative part of the advection is added in updateConservative.
     */
    computeFluxes() {
        const U = this.U;
        const nx = this.nx;
        const nvar = this.nvar;
        const M = this.numMaterials;
        const p = this.pc;
        const rc2 = this.rc2;
        const tiny = 1e-10;
        let maxWaveSpeed = 0;

        if (this.interfaceMethod === 'mixed') {
            this.updateEOS(U);
        }
        if (!this.tabulated) {
            // Pressure and rho*a^2 = gamma*(p + pinf) = (p*(G+1) + P)/G
            const G = this.eosG, P = this.eosP;
            for (let i = 0; i < nx; i++) {
                const c = i * nvar;
                const rho = Math.max(U[c], tiny);
                const u = U[c + 1] / rho;
                p[i] = (U[c + 2] - 0.5 * rho * u * u - P[i]) / G[i];
                rc2[i] = (p[i] * (G[i] + 1) + P[i]) / G[i];
                if (!(rc2[i] > tiny)) {
                    rc2[i] = tiny;
                    p[i] = (tiny * G[i] - P[i]) / (G[i] + 1);
                }
            }
        }

        for (let i = 0; i <= nx; i++) {
            const f = i * nvar;
            let waveSpeed;

            if (i === 0) {
                // Left wall: mirror state with reversed velocity
                waveSpeed = this.hllcFlux(U[0], -U[1], U[2], p[0], rc2[0],
                                          U[0], U[1], U[2], p[0], rc2[0], f);
            } else if (i === nx) {
                // Right wall
                const c = (nx - 1) * nvar;
                waveSpeed = this.hllcFlux(U[c], U[c + 1], U[c + 2], p[nx - 1], rc2[nx - 1],
                                          U[c], -U[c + 1], U[c + 2], p[nx - 1], rc2[nx - 1], f);
            } else {
                const l = (i - 1) * nvar;
                const r = i * nvar;
                waveSpeed = this.hllcFlux(U[l], U[l + 1], U[l + 2], p[i - 1], rc2[i - 1],
                                          U[r], U[r + 1], U[r + 2], p[i], rc2[i], f);
            }

            const uFace = this._uFace;
            this.faceVelocity[i] = uFace;

            if (M > 0) {
                // Upwind volume fractions (walls: the single interior cell)
                const cell = i === 0 ? 0 : (i === nx ? nx - 1 : (this._upwind === 0 ? i - 1 : i));
                const src = cell * nvar + 3;
                for (let m = 0; m < M; m++) {
                    this.flux[f + 3 + m] = U[src + m] * uFace;
                }
                if (this.tabulated) {
                    // Partial densities ride on the mass flux with the upwind mass fraction
                    const massFlux = this.flux[f];
                    const rhoUp = Math.max(U[cell * nvar], tiny);
                    for (let m = 0; m < M; m++) {
                        this.flux[f + 3 + M + m] = massFlux * U[src + M + m] / rhoUp;
                    }
                }
            }

            if (waveSpeed > maxWaveSpeed) maxWaveSpeed = waveSpeed;
        }

        return maxWaveSpeed;
    }

    /**
     * Forward-Euler update Uout = U - dt * dF/dx using the stored fluxes.
     * Volume fractions obey the advection equation d(alpha)/dt + u d(alpha)/dx = 0,
     * discretized as div(alpha*u) - alpha*div(u) with the same face velocities
     * as the HLLC flux so that material interfaces stay free of pressure
     * oscillations.
     */
    updateConservative(U, Uout, dt) {
        const nvar = this.nvar;
        const lambda = dt / this.dx;
        const flux = this.flux;

        for (let i = 0; i < this.nx; i++) {
            const idx = i * nvar;
            const fL = idx;
            const fR = idx + nvar;

            for (let k = 0; k < 3; k++) {
                Uout[idx + k] = U[idx + k] - lambda * (flux[fR + k] - flux[fL + k]);
            }

            if (nvar > 3) {
                const lastAlpha = 3 + this.numMaterials;
                const divU = this.faceVelocity[i + 1] - this.faceVelocity[i];
                for (let k = 3; k < lastAlpha; k++) {
                    const alpha = U[idx + k];
                    Uout[idx + k] = alpha - lambda * (flux[fR + k] - flux[fL + k] - alpha * divU);
                }
                // Partial densities are conserved
                for (let k = lastAlpha; k < nvar; k++) {
                    Uout[idx + k] = U[idx + k] - lambda * (flux[fR + k] - flux[fL + k]);
                }
            }
        }
    }

    /**
     * Update per-cell gas properties after a completed time step
     */
    updateGasProperties() {
        if (this.interfaceMethod === 'mixed') {
            // Remove round-off drift: clamp to [0, 1] and renormalize
            const nvar = this.nvar;
            const M = this.numMaterials;
            for (let i = 0; i < this.nx; i++) {
                const base = i * nvar + 3;
                let total = 0;
                for (let m = 0; m < M; m++) {
                    const alpha = Math.min(1, Math.max(0, this.U[base + m]));
                    this.U[base + m] = alpha;
                    total += alpha;
                }
                if (total > 1e-12) {
                    for (let m = 0; m < M; m++) this.U[base + m] /= total;
                } else {
                    for (let m = 0; m < M; m++) this.U[base + m] = 1 / M;
                }
                if (this.tabulated) {
                    for (let m = 0; m < M; m++) {
                        if (this.U[base + M + m] < 0) this.U[base + M + m] = 0;
                    }
                }
            }
        }
        this.updateEOS(this.U);
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
            
            // Store trajectory point with velocity
            tracer.trajectory.push({t: this.t, x: tracer.x, v: velocity});
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
            pinf: new Float64Array(this.pinf),  // Stiffening pressure (0 for ideal gas)
            mw: new Float64Array(this.mw),      // Molecular weight
            gasId: Array.from(this.gasId)       // Gas identifier
        };
        
        // Include material fractions for mixed-cell method
        if (this.interfaceMethod === 'mixed' && this.materialFractions) {
            snapshot.materialFractions = new Float64Array(this.materialFractions);
            snapshot.numMaterials = this.numMaterials;
            snapshot.materialProperties = this.materialProperties.map(prop => ({
                gamma: prop.gamma, mw: prop.mw, pinf: prop.pinf, gasId: prop.gasId
            }));
        }
        
        this.xtData.push(snapshot);
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
            timeSteps: this.timeSteps,
            wallTime: this.wallTime || 0,
            interfaceMethod: this.interfaceMethod
        };
    }
    
    /**
     * Advance the simulation for at most budgetMs of wall time.
     * @returns {boolean} true once the final time has been reached
     */
    advance(budgetMs) {
        const deadline = Date.now() + budgetMs;
        const begin = Date.now();
        try {
            let sinceClock = 0;
            while (this.t < this.finalTime) {
                this.integrator.step(this);

                if (!isFinite(this.dt) || this.dt <= 0) {
                    throw new Error(`Invalid time step detected: dt = ${this.dt}`);
                }
                if (this.timeSteps > 1000000) {
                    throw new Error('Exceeded 1,000,000 time steps. Reduce the simulated time or the grid resolution.');
                }
                // Reading the clock every step would dominate small runs
                if (++sinceClock >= 20) {
                    sinceClock = 0;
                    if (Date.now() >= deadline) break;
                }
            }
        } catch (error) {
            throw new Error(`Euler solver failed at t = ${(this.t * 1000).toFixed(4)} ms after ${this.timeSteps} steps: ${error.message}`);
        } finally {
            this.wallTime = (this.wallTime || 0) + (Date.now() - begin) / 1000;
        }
        return this.t >= this.finalTime;
    }

    /**
     * Run simulation to completion (blocking)
     */
    run(progressCallback) {
        while (!this.advance(100)) {
            if (progressCallback) progressCallback(this.t / this.finalTime);
        }
        if (progressCallback) progressCallback(1.0);
    }
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { EulerSolver };
}
