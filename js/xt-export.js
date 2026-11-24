/**
 * HDF5 export of x-t simulation results.
 * Relies on globals from app-xt.js: xtResults, xtRunConfig, showXTMessage.
 */

const H5WASM_URL = 'https://cdn.jsdelivr.net/npm/h5wasm@0.7.3/dist/esm/hdf5_hl.js';
let h5wasmPromise = null;

function loadH5Wasm() {
    if (!h5wasmPromise) {
        h5wasmPromise = import(H5WASM_URL).then(function (mod) {
            return mod.default.ready.then(function () { return mod.default; });
        }).catch(function (err) {
            h5wasmPromise = null;   // allow a retry
            throw new Error('could not load the HDF5 library (' + err.message + ')');
        });
    }
    return h5wasmPromise;
}

async function exportXTDataHDF5() {
    if (!xtResults) {
        showXTMessage('No simulation results to export.', 'warning');
        return;
    }
    
    try {
        // h5wasm (~1 MB) is only fetched the first time an HDF5 export is requested
        showXTMessage('Loading HDF5 library…', 'info');
        const h5wasm = await loadH5Wasm();
        const { FS } = await h5wasm.ready;
        
        // Generate datetime-stamped filename
        const now = new Date();
        const year = now.getFullYear();
        const month = String(now.getMonth() + 1).padStart(2, '0');
        const day = String(now.getDate()).padStart(2, '0');
        const hour = String(now.getHours()).padStart(2, '0');
        const minute = String(now.getMinutes()).padStart(2, '0');
        const second = String(now.getSeconds()).padStart(2, '0');
        const filename = `xtcalc_${year}${month}${day}_${hour}${minute}${second}.h5`;
        
        showXTMessage('Creating HDF5 file...', 'info');
        
        // Create HDF5 file
        const f = new h5wasm.File(filename, 'w');
        
        // Extract dimensions from full data
        const num_snapshots = xtResults.xtData.length;
        const nx = xtResults.x.length;
        
        // Check if we have mixed-cell data to determine species
        const hasMixedCellData = xtResults.xtData[0].materialFractions !== undefined;
        
        let speciesArray;
        let n_species;
        
        if (hasMixedCellData) {
            // For mixed-cell: species are the original materials, not computed mixtures
            const firstSnapshot = xtResults.xtData[0];
            speciesArray = firstSnapshot.materialProperties.map(prop => ({
                gamma: prop.gamma,
                mw: prop.mw,
                pinf: prop.pinf || 0,
                names: [prop.gasId],
                gasIds: [prop.gasId]
            }));
            n_species = speciesArray.length;
        } else {
            // For sharp interface: deduplicate species by (gamma, MW) from cell data
            const speciesMap = new Map(); // key: "gamma_mw" -> {gamma, mw, names: [], gasIds: []}
            
            // Collect all unique (gamma, MW) pairs from full data
            for (let i = 0; i < num_snapshots; i++) {
                const snapshot = xtResults.xtData[i];
                for (let j = 0; j < nx; j++) {
                    const gamma = snapshot.gamma[j];
                    const mw = snapshot.mw[j];
                    const pinf = snapshot.pinf ? snapshot.pinf[j] : 0;
                    const gasId = snapshot.gasId[j];
                    
                    // Create key with tolerance
                    const key = `${gamma.toFixed(6)}_${mw.toFixed(4)}`;
                    
                    if (!speciesMap.has(key)) {
                        speciesMap.set(key, {
                            gamma: gamma,
                            mw: mw,
                            pinf: pinf,
                            names: [gasId],
                            gasIds: [gasId]
                        });
                    } else {
                        const species = speciesMap.get(key);
                        if (!species.gasIds.includes(gasId)) {
                            species.names.push(gasId);
                            species.gasIds.push(gasId);
                        }
                    }
                }
            }
            
            speciesArray = Array.from(speciesMap.values());
            n_species = speciesArray.length;
        }
        
        // Create groups
        f.create_group('metadata');
        f.create_group('species');
        f.create_group('coordinates');
        f.create_group('fields');
        f.create_group('tracers');
        
        // Write metadata
        f.get('metadata').create_dataset({name: 'nx', data: [nx], dtype: '<i'});
        f.get('metadata').create_dataset({name: 'num_snapshots', data: [num_snapshots], dtype: '<i'});
        f.get('metadata').create_dataset({name: 'n_species', data: [n_species], dtype: '<i'});
        f.get('metadata').create_dataset({name: 'cfl', data: [xtRunConfig.cfl], dtype: '<f'});
        f.get('metadata').create_dataset({name: 'final_time', data: [xtResults.t], dtype: '<d'});
        f.get('metadata').create_dataset({name: 'integrator', data: [xtRunConfig.integrator]});
        f.get('metadata').create_dataset({name: 'interface_method', data: [xtResults.interfaceMethod]});
        f.get('metadata').create_dataset({name: 'gas_model', data: [xtRunConfig.gasModel]});
        f.get('metadata').create_dataset({name: 'spatial_order', data: [xtRunConfig.order], dtype: '<i'});
        f.get('metadata').create_dataset({name: 'simulation_datetime', data: [now.toISOString()]});
        f.get('metadata').create_dataset({name: 'wall_time', data: [xtResults.wallTime || 0], dtype: '<d'});
        
        // Write species information
        const gammaArray = speciesArray.map(s => s.gamma);
        const mwArray = speciesArray.map(s => s.mw);
        const namesArray = speciesArray.map(s => s.names.join(' / '));
        const gasIdsArray = speciesArray.map(s => s.gasIds[0]);
        
        f.get('species').create_dataset({name: 'gamma', data: gammaArray, dtype: '<d'});
        f.get('species').create_dataset({name: 'molecular_weight', data: mwArray, dtype: '<d'});
        // Stiffened-gas pressure constant; zero for ideal gases
        f.get('species').create_dataset({name: 'stiffening_pressure', data: speciesArray.map(s => s.pinf || 0), dtype: '<d'});
        f.get('species').create_dataset({name: 'names', data: namesArray});
        f.get('species').create_dataset({name: 'representative_gas_ids', data: gasIdsArray});
        
        // Write coordinates
        f.get('coordinates').create_dataset({name: 'x', data: xtResults.x, dtype: '<d'});
        const timeArray = xtResults.xtData.map(snap => snap.t);
        f.get('coordinates').create_dataset({name: 't', data: timeArray, dtype: '<d'});
        
        // Build 2D arrays for fields
        const density2D = new Float64Array(num_snapshots * nx);
        const velocity2D = new Float64Array(num_snapshots * nx);
        const pressure2D = new Float64Array(num_snapshots * nx);
        const temperature2D = new Float64Array(num_snapshots * nx);
        const gamma2D = new Float64Array(num_snapshots * nx);
        const mw2D = new Float64Array(num_snapshots * nx);
        
        // Build 3D volume fraction array (num_snapshots x nx x n_species)
        const volumeFraction3D = new Float32Array(num_snapshots * nx * n_species);
        volumeFraction3D.fill(0);
        
        for (let i = 0; i < num_snapshots; i++) {
            const snapshot = xtResults.xtData[i];
            for (let j = 0; j < nx; j++) {
                const idx2D = i * nx + j;
                density2D[idx2D] = snapshot.rho[j];
                velocity2D[idx2D] = snapshot.u[j];
                pressure2D[idx2D] = snapshot.p[j];
                temperature2D[idx2D] = snapshot.T[j];
                gamma2D[idx2D] = snapshot.gamma[j];
                mw2D[idx2D] = snapshot.mw[j];
                
                if (hasMixedCellData) {
                    // Mixed-cell method: use actual material fractions
                    const numMaterials = snapshot.numMaterials;
                    
                    // Map each material to its corresponding species index
                    for (let m = 0; m < numMaterials; m++) {
                        const matProps = snapshot.materialProperties[m];
                        const matKey = `${matProps.gamma.toFixed(6)}_${matProps.mw.toFixed(4)}`;
                        
                        // Find matching species
                        let speciesIdx = 0;
                        for (let k = 0; k < speciesArray.length; k++) {
                            const speciesKey = `${speciesArray[k].gamma.toFixed(6)}_${speciesArray[k].mw.toFixed(4)}`;
                            if (speciesKey === matKey) {
                                speciesIdx = k;
                                break;
                            }
                        }
                        
                        // Get material fraction for this cell and material
                        const matFracIdx = j * numMaterials + m;
                        const alpha = snapshot.materialFractions[matFracIdx];
                        
                        // Store in 3D array
                        const idx3D = i * nx * n_species + j * n_species + speciesIdx;
                        volumeFraction3D[idx3D] = alpha;
                    }
                } else {
                    // Sharp method: single species per cell
                    const cellGamma = snapshot.gamma[j];
                    const cellMw = snapshot.mw[j];
                    const key = `${cellGamma.toFixed(6)}_${cellMw.toFixed(4)}`;
                    
                    let speciesIdx = 0;
                    for (let k = 0; k < speciesArray.length; k++) {
                        const speciesKey = `${speciesArray[k].gamma.toFixed(6)}_${speciesArray[k].mw.toFixed(4)}`;
                        if (speciesKey === key) {
                            speciesIdx = k;
                            break;
                        }
                    }
                    
                    // Set volume fraction to 1.0 for this species
                    const idx3D = i * nx * n_species + j * n_species + speciesIdx;
                    volumeFraction3D[idx3D] = 1.0;
                }
            }
        }
        
        // Write 2D field datasets with compression
        const chunkSize2D = [Math.min(100, num_snapshots), Math.min(100, nx)];
        
        f.get('fields').create_dataset({
            name: 'density', 
            data: density2D, 
            shape: [num_snapshots, nx], 
            dtype: '<d',
            chunks: chunkSize2D,
            compression: 6
        });
        
        f.get('fields').create_dataset({
            name: 'velocity', 
            data: velocity2D, 
            shape: [num_snapshots, nx], 
            dtype: '<d',
            chunks: chunkSize2D,
            compression: 6
        });
        
        f.get('fields').create_dataset({
            name: 'pressure', 
            data: pressure2D, 
            shape: [num_snapshots, nx], 
            dtype: '<d',
            chunks: chunkSize2D,
            compression: 6
        });
        
        f.get('fields').create_dataset({
            name: 'temperature', 
            data: temperature2D, 
            shape: [num_snapshots, nx], 
            dtype: '<d',
            chunks: chunkSize2D,
            compression: 6
        });
        
        f.get('fields').create_dataset({
            name: 'gamma', 
            data: gamma2D, 
            shape: [num_snapshots, nx], 
            dtype: '<d',
            chunks: chunkSize2D,
            compression: 6
        });
        
        f.get('fields').create_dataset({
            name: 'molecular_weight', 
            data: mw2D, 
            shape: [num_snapshots, nx], 
            dtype: '<d',
            chunks: chunkSize2D,
            compression: 6
        });
        
        // Write 3D volume fraction dataset with compression
        const chunkSize3D = [Math.min(100, num_snapshots), Math.min(100, nx), n_species];
        
        f.get('fields').create_dataset({
            name: 'volume_fraction', 
            data: volumeFraction3D, 
            shape: [num_snapshots, nx, n_species], 
            dtype: '<f',
            chunks: chunkSize3D,
            compression: 6
        });
        
        // Write tracer data
        f.get('tracers').create_dataset({name: 'num_tracers', data: [xtResults.tracers.length], dtype: '<i'});
        
        for (let i = 0; i < xtResults.tracers.length; i++) {
            const tracer = xtResults.tracers[i];
            const groupName = 'interface_' + i;
            f.get('tracers').create_group(groupName);
            
            const times = tracer.trajectory.map(p => p.t);
            const positions = tracer.trajectory.map(p => p.x);
            const velocities = tracer.trajectory.map(p => p.v);
            
            f.get('tracers/' + groupName).create_dataset({name: 'times', data: times, dtype: '<d'});
            f.get('tracers/' + groupName).create_dataset({name: 'positions', data: positions, dtype: '<d'});
            f.get('tracers/' + groupName).create_dataset({name: 'velocities', data: velocities, dtype: '<d'});
        }
        
        // Flush and close file
        f.flush();
        f.close();
        
        // Read file from Emscripten FS and download
        const fileData = FS.readFile(filename);
        const blob = new Blob([fileData], { type: 'application/x-hdf5' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.click();
        URL.revokeObjectURL(url);
        
        // Clean up
        FS.unlink(filename);
        
        // Show file size info
        const sizeKB = (fileData.length / 1024).toFixed(1);
        showXTMessage(`HDF5 file exported successfully! (${sizeKB} KB, compressed)`, 'success');
        
    } catch (error) {
        showXTMessage('Error exporting HDF5: ' + error.message, 'error');
        console.error('HDF5 export error:', error);
    }
}

