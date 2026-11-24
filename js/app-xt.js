// X-T Diagram state
var xtSlabs = [];
var xtSlabCounter = 0;
var xtSlabMixtures = {}; // Store mixture compositions for slabs with custom mixtures
var xtSolver = null;
var xtVisualization = null;
var xtResults = null;
var xtSuppressRender = false; // Flag to suppress rendering during batch operations

// Update interface method information based on selection
function updateInterfaceMethodInfo() {
    var method = document.getElementById('xtInterfaceMethod').value;
    var infoSpan = document.getElementById('interfaceMethodInfo');
    
    if (method === 'sharp') {
        infoSpan.innerHTML = '<span style="color:#6c757d;">Fastest, may show oscillations near interfaces</span>';
    } else if (method === 'mixed') {
        infoSpan.innerHTML = '<span style="color:#dc3545;">Most accurate, eliminates oscillations (~15% slower)</span>';
    }
}

// Toggle X-T Diagram section
function toggleXTDiagramSection() {
    var section = document.getElementById('xtDiagramSection');
    var btn = document.getElementById('toggleXTDiagram');
    
    if (section.style.display === 'none') {
        section.style.display = 'block';
        btn.textContent = 'Hide X-T Diagram Calculator ▲';
        
        // Initialize visualization if not already done
        if (!xtVisualization) {
            xtVisualization = new XTVisualization('xtCanvas', {
                width: 1000,
                height: 600
            });
        }
        
        // Add default slabs if none exist
        if (xtSlabs.length === 0) {
            addXTSlab({ gas: '4', pressure: 400000, temperature: 300, length: 3 }); // Driver (Helium)
            addXTSlab({ gas: '1', pressure: 101325, temperature: 300, length: 6 }); // Driven (Air)
        }
    } else {
        section.style.display = 'none';
        btn.textContent = 'Show X-T Diagram Calculator ▼';
    }
}

// Use current calculation to populate slabs
function useCurrentCalculation() {
    // Check if calculation has been run
    var pfour = document.Shock.pfour.value;
    if (!pfour || pfour === '') {
        document.getElementById('xtImportWarning').textContent = '⚠ Please run the driver/driven calculation first.';
        document.getElementById('xtImportInfo').textContent = '';
        return;
    }
    
    // Clear warnings and info
    document.getElementById('xtImportWarning').textContent = '';
    document.getElementById('xtImportInfo').textContent = '';
    
// Preserve existing slab lengths for slabs #1 and #2
var existingDriverLength = xtSlabs[0] ? xtSlabs[0].length : 3; // Default 3m for driver
var existingDrivenLength = xtSlabs[1] ? xtSlabs[1].length : 6; // Default 6m for driven

// Store any additional slabs beyond #2
var additionalSlabs = xtSlabs.slice(2);

// Clear only the first two slabs
xtSlabs = [];
xtSlabCounter = 0;
xtSlabMixtures = {}; // Clear mixture data

    // Suppress rendering during batch import
    xtSuppressRender = true;
    
    // Get driver gas info
    var driverGasId = document.Shock.drivergas.value;
    var driverP = parseFloat(pfour) * 1e6; // Convert MPa to Pa
    var driverT = parseFloat(document.Shock.tid.value);
    
    // Get driven gas info
    var drivenGasId = document.Shock.drivengas.value;
    var drivenP = parseFloat(document.Shock.pinit.value);
    var drivenT = parseFloat(document.Shock.tinit.value);
    
    // Add driver slab with preserved length
    if (driverGasId === 'custom') {
        // Import driver custom mixture - capture ID AFTER adding slab
        addXTSlab({ gas: 'custom_mixture', pressure: driverP, temperature: driverT, length: existingDriverLength });
        var driverSlabId = xtSlabs[xtSlabs.length - 1].id;
        
        // Copy mixture composition from shock tube calculator
        xtSlabMixtures[driverSlabId] = [];
        for (var i = 0; i < mixtures.driver.length; i++) {
            var index = mixtures.driver[i];
            var gasSelect = document.getElementById('driver_gas_' + index);
            var fracInput = document.getElementById('driver_frac_' + index);
            
            if (gasSelect && fracInput && gasSelect.value) {
                var compId = xtMixtureComponentCounter++;
                xtSlabMixtures[driverSlabId].push({
                    id: compId,
                    gasId: gasSelect.value,
                    fraction: parseFloat(fracInput.value)
                });
            }
        }
    } else {
        addXTSlab({ gas: driverGasId, pressure: driverP, temperature: driverT, length: existingDriverLength });
    }
    
    // Add driven slab with preserved length
    if (drivenGasId === 'custom') {
        // Import driven custom mixture - capture ID AFTER adding slab
        addXTSlab({ gas: 'custom_mixture', pressure: drivenP, temperature: drivenT, length: existingDrivenLength });
        var drivenSlabId = xtSlabs[xtSlabs.length - 1].id;
        
        // Copy mixture composition from shock tube calculator
        xtSlabMixtures[drivenSlabId] = [];
        for (var i = 0; i < mixtures.driven.length; i++) {
            var index = mixtures.driven[i];
            var gasSelect = document.getElementById('driven_gas_' + index);
            var fracInput = document.getElementById('driven_frac_' + index);
            
            if (gasSelect && fracInput && gasSelect.value) {
                var compId = xtMixtureComponentCounter++;
                xtSlabMixtures[drivenSlabId].push({
                    id: compId,
                    gasId: gasSelect.value,
                    fraction: parseFloat(fracInput.value)
                });
            }
        }
    } else {
        addXTSlab({ gas: drivenGasId, pressure: drivenP, temperature: drivenT, length: existingDrivenLength });
    }
    
    // Restore any additional slabs beyond #2
    for (var i = 0; i < additionalSlabs.length; i++) {
        xtSlabs.push(additionalSlabs[i]);
    }
    
    // Re-enable rendering and render once with all data ready
    xtSuppressRender = false;
    renderXTSlabs();
    
    // Show import info
    var driverDesc = driverGasId === 'custom' ? 'Custom Mixture' : getAllGases()[driverGasId].name;
    var drivenDesc = drivenGasId === 'custom' ? 'Custom Mixture' : getAllGases()[drivenGasId].name;
    var info = '✓ Imported: ' + driverDesc + ' (driver) + ' + drivenDesc + ' (driven)';
    document.getElementById('xtImportInfo').textContent = info;
}

// Add a gas slab
function addXTSlab(config = {}) {
    var id = xtSlabCounter++;
    
    // Get defaults from preceding slab if it exists
    var defaultPressure = 101325;  // 1 atm
    var defaultTemperature = 300;   // K
    
    if (xtSlabs.length > 0) {
        var precedingSlab = xtSlabs[xtSlabs.length - 1];
        defaultPressure = precedingSlab.pressure;
        defaultTemperature = precedingSlab.temperature;
    }
    
    var slab = {
        id: id,
        gas: config.gas || '',
        pressure: config.pressure || defaultPressure,
        temperature: config.temperature || defaultTemperature,
        length: config.length || 1
    };
    
    xtSlabs.push(slab);
    
    // Only render if not suppressed (during batch operations)
    if (!xtSuppressRender) {
        renderXTSlabs();
    }
}

// Remove a slab
function removeXTSlab(id) {
    xtSlabs = xtSlabs.filter(s => s.id !== id);
    renderXTSlabs();
}

// Move slab up
function moveXTSlabUp(id) {
    var index = xtSlabs.findIndex(s => s.id === id);
    if (index > 0) {
        var temp = xtSlabs[index - 1];
        xtSlabs[index - 1] = xtSlabs[index];
        xtSlabs[index] = temp;
        renderXTSlabs();
    }
}

// Move slab down
function moveXTSlabDown(id) {
    var index = xtSlabs.findIndex(s => s.id === id);
    if (index < xtSlabs.length - 1) {
        var temp = xtSlabs[index + 1];
        xtSlabs[index + 1] = xtSlabs[index];
        xtSlabs[index] = temp;
        renderXTSlabs();
    }
}

// Update slab property
function updateXTSlab(id, property, value) {
    var slab = xtSlabs.find(s => s.id === id);
    if (slab) {
        slab[property] = parseFloat(value) || value;
    }
}

// Generate gas options HTML
function generateXTGasOptions(selectedGas) {
    var allGases = getAllGases();
    var options = '<option value="">Select Gas</option>';
    options += '<option value="custom_mixture"' + (selectedGas === 'custom_mixture' ? ' selected' : '') + '>Custom Mixture</option>';
    options += '<option disabled>──────────</option>';
    
    var gasEntries = Object.entries(allGases).sort(function(a, b) {
        return a[1].name.localeCompare(b[1].name);
    });
    
    for (var i = 0; i < gasEntries.length; i++) {
        var id = gasEntries[i][0];
        var gas = gasEntries[i][1];
        
        // Skip Acetone Vapor in single gas selection (only available in custom mixtures)
        if (gas.name === 'Acetone Vapor') continue;
        
        var selected = id === selectedGas ? 'selected' : '';
        var prefix = id.startsWith('custom_') ? '★ ' : '';
        options += '<option value="' + id + '" ' + selected + '>' + prefix + gas.name + '</option>';
    }
    
    return options;
}

// Render slab list
function renderXTSlabs() {
    var container = document.getElementById('xtSlabList');
    
    if (xtSlabs.length === 0) {
        container.innerHTML = '<p style="color: #6c757d;">No slabs configured. Add at least one slab to begin.</p>';
        return;
    }
    
    // Store current values before rebuilding
    var currentValues = {};
    for (var i = 0; i < xtSlabs.length; i++) {
        var slab = xtSlabs[i];
        var selectElem = document.getElementById('xt_slab_gas_' + slab.id);
        if (selectElem) {
            currentValues[slab.id] = selectElem.value;
        }
    }
    
    var html = '';
    for (var i = 0; i < xtSlabs.length; i++) {
        var slab = xtSlabs[i];
        // Use stored value if available, otherwise use slab.gas
        var selectedGas = currentValues[slab.id] !== undefined ? currentValues[slab.id] : slab.gas;
        
        html += '<div class="slab-item" id="xt_slab_item_' + slab.id + '">';
        html += '<div class="slab-number">#' + (i + 1) + '</div>';
        html += '<div class="slab-controls">';
        html += '<div class="slab-field">';
        html += '<label>Gas</label>';
        html += '<select id="xt_slab_gas_' + slab.id + '" onchange="updateXTSlab(' + slab.id + ', \'gas\', this.value); toggleXTSlabMixture(' + slab.id + ')">';
        html += generateXTGasOptions(selectedGas);
        html += '</select>';
        html += '</div>';
        html += '<div class="slab-field">';
        html += '<label>Pressure (kPa)</label>';
        html += '<input type="number" value="' + (slab.pressure / 1000).toFixed(1) + '" ';
        html += 'onchange="updateXTSlab(' + slab.id + ', \'pressure\', this.value * 1000)" step="1" min="1">';
        html += '</div>';
        html += '<div class="slab-field">';
        html += '<label>Temperature (K)</label>';
        html += '<input type="number" value="' + slab.temperature + '" ';
        html += 'onchange="updateXTSlab(' + slab.id + ', \'temperature\', this.value)" step="1" min="1">';
        html += '</div>';
        html += '<div class="slab-field">';
        html += '<label>Length (m)</label>';
        html += '<input type="number" value="' + slab.length + '" ';
        html += 'onchange="updateXTSlab(' + slab.id + ', \'length\', this.value)" step="0.1" min="0.1">';
        html += '</div>';
        html += '</div>';
        html += '<div class="slab-actions">';
        html += '<button type="button" onclick="moveXTSlabUp(' + slab.id + ')" ' + (i === 0 ? 'disabled' : '') + '>↑</button>';
        html += '<button type="button" onclick="moveXTSlabDown(' + slab.id + ')" ' + (i === xtSlabs.length - 1 ? 'disabled' : '') + '>↓</button>';
        html += '<button type="button" class="danger" onclick="removeXTSlab(' + slab.id + ')">Remove</button>';
        html += '</div>';
        html += '</div>';
    }
    
    container.innerHTML = html;
    
    // Explicitly set dropdown values after HTML render (more reliable than selected attribute)
    for (var i = 0; i < xtSlabs.length; i++) {
        var slab = xtSlabs[i];
        var select = document.getElementById('xt_slab_gas_' + slab.id);
        if (select) {
            select.value = slab.gas;
        }
    }
    
    // Re-initialize mixture UIs for slabs that have custom mixture selected
    for (var i = 0; i < xtSlabs.length; i++) {
        var slab = xtSlabs[i];
        if (slab.gas === 'custom_mixture') {
            toggleXTSlabMixture(slab.id);
        }
    }
}

// Validate configuration
function validateXTConfiguration() {
    if (xtSlabs.length === 0) {
        showXTMessage('Please add at least one gas slab.', 'error');
        return false;
    }
    
    // Check for same-name gases with different properties
    var gasNameMap = {}; // name → {gamma, mw, slabIndex}
    
    for (var i = 0; i < xtSlabs.length; i++) {
        var slab = xtSlabs[i];
        
        if (!slab.gas) {
            showXTMessage('Slab #' + (i + 1) + ': Please select a gas.', 'error');
            return false;
        }
        
        if (slab.pressure <= 0) {
            showXTMessage('Slab #' + (i + 1) + ': Pressure must be positive.', 'error');
            return false;
        }
        
        if (slab.temperature <= 0) {
            showXTMessage('Slab #' + (i + 1) + ': Temperature must be positive.', 'error');
            return false;
        }
        
        if (slab.length <= 0) {
            showXTMessage('Slab #' + (i + 1) + ': Length must be positive.', 'error');
            return false;
        }
        
        // Validate gas properties for conflicts
        try {
            var gasProps = getXTSlabGasProperties(slab);
            var name = gasProps.gasName;
            
            if (gasNameMap[name]) {
                // Name exists, check if properties match
                var existing = gasNameMap[name];
                var gammaMatch = Math.abs(existing.gamma - gasProps.gamma) < 1e-6;
                var mwMatch = Math.abs(existing.mw - gasProps.mw) < 1e-6;
                
                if (!gammaMatch || !mwMatch) {
                    var errorMsg = 'Conflicting gas properties detected for "' + name + '":\n\n' +
                                   'Slab #' + (i + 1) + ' has: γ=' + gasProps.gamma.toFixed(6) + ', MW=' + gasProps.mw.toFixed(4) + ' g/mol\n' +
                                   'Slab #' + existing.slabIndex + ' has: γ=' + existing.gamma.toFixed(6) + ', MW=' + existing.mw.toFixed(4) + ' g/mol\n\n' +
                                   'Please ensure all gases with the same name have identical thermodynamic properties, or rename one of them to avoid confusion.';
                    showXTMessage(errorMsg, 'error');
                    return false;
                }
            } else {
                gasNameMap[name] = {
                    gamma: gasProps.gamma,
                    mw: gasProps.mw,
                    slabIndex: i + 1
                };
            }
        } catch (error) {
            showXTMessage('Slab #' + (i + 1) + ': ' + error.message, 'error');
            return false;
        }
    }
    
    return true;
}

// Run X-T simulation
function runXTSimulation() {
    if (!validateXTConfiguration()) {
        return;
    }
    
    var runBtn = document.getElementById('xtRunBtn');
    runBtn.disabled = true;
    runBtn.textContent = 'Running...';
    
    document.getElementById('xtProgressContainer').style.display = 'block';
    document.getElementById('xtCanvasContainer').style.display = 'none';
    document.getElementById('xtStatsContainer').style.display = 'none';
    
    var integrator = document.getElementById('xtIntegrator').value;
    var interfaceMethod = document.getElementById('xtInterfaceMethod').value;
    var nx = parseInt(document.getElementById('xtGridPoints').value);
    var simTime = parseFloat(document.getElementById('xtSimTime').value) / 1000;
    var snapshotInterval = parseFloat(document.getElementById('xtSnapshotInterval').value) / 1000;
    var cfl = parseFloat(document.getElementById('xtCflNumber').value);
    
    var slabData = [];
    try {
        slabData = xtSlabs.map(function(slab) {
            var props = getXTSlabGasProperties(slab);
            return {
                gasId: props.gasName,  // Use human-readable name
                gamma: props.gamma,
                mw: props.mw,
                pressure: slab.pressure,
                temperature: slab.temperature,
                length: slab.length
            };
        });
    } catch (error) {
        showXTMessage('Configuration error: ' + error.message, 'error');
        runBtn.disabled = false;
        runBtn.textContent = 'Run Simulation';
        return;
    }
    
    xtSolver = new EulerSolver({
        nx: nx,
        finalTime: simTime,
        cfl: cfl,
        xtInterval: snapshotInterval,
        integrator: integrator,
        interfaceMethod: interfaceMethod
    });
    
    try {
        xtSolver.initialize(slabData);
        showXTMessage('Simulation initialized. Running solver...', 'info');
        
        setTimeout(function() {
            try {
                // Capture wall time
                var wallStartTime = performance.now();
                
                xtSolver.run(function(progress) {
                    var percent = Math.round(progress * 100);
                    document.getElementById('xtProgressFill').style.width = percent + '%';
                    document.getElementById('xtProgressFill').textContent = percent + '%';
                });
                
                var wallEndTime = performance.now();
                xtSolver.wallTime = (wallEndTime - wallStartTime) / 1000; // Convert to seconds
                
                xtResults = xtSolver.getResults();
                xtVisualization.setData(xtResults);
                xtVisualization.render();
                
                document.getElementById('xtCanvasContainer').style.display = 'block';
                displayXTStatistics(xtResults);
                
                showXTMessage('Simulation completed successfully!', 'success');
            } catch (error) {
                showXTMessage('Simulation error: ' + error.message, 'error');
                console.error(error);
            } finally {
                runBtn.disabled = false;
                runBtn.textContent = 'Run Simulation';
                document.getElementById('xtProgressContainer').style.display = 'none';
            }
        }, 100);
    } catch (error) {
        showXTMessage('Initialization error: ' + error.message, 'error');
        console.error(error);
        runBtn.disabled = false;
        runBtn.textContent = 'Run Simulation';
        document.getElementById('xtProgressContainer').style.display = 'none';
    }
}

// Display statistics
function displayXTStatistics(results) {
    var container = document.getElementById('xtStatsContent');
    var totalLength = results.x[results.x.length - 1];
    var numSnapshots = results.xtData.length;
    var wallTime = results.wallTime || 0;
    
    container.innerHTML = 
        '<div><div style="font-weight:600;color:#495057;">Time Steps</div><div style="font-size:18px;">' + results.timeSteps.toLocaleString() + '</div></div>' +
        '<div><div style="font-weight:600;color:#495057;">Final Time</div><div style="font-size:18px;">' + (results.t * 1000).toFixed(2) + ' ms</div></div>' +
        '<div><div style="font-weight:600;color:#495057;">Wall Time</div><div style="font-size:18px;">' + wallTime.toFixed(2) + ' s</div></div>' +
        '<div><div style="font-weight:600;color:#495057;">X-T Snapshots</div><div style="font-size:18px;">' + numSnapshots + '</div></div>' +
        '<div><div style="font-weight:600;color:#495057;">Domain Length</div><div style="font-size:18px;">' + totalLength.toFixed(2) + ' m</div></div>';
    
    document.getElementById('xtStatsContainer').style.display = 'block';
}

// Reset simulation
function resetXTSimulation() {
    xtSolver = null;
    xtResults = null;
    document.getElementById('xtCanvasContainer').style.display = 'none';
    document.getElementById('xtStatsContainer').style.display = 'none';
    document.getElementById('xtMessageContainer').innerHTML = '';
    document.getElementById('xtProgressContainer').style.display = 'none';
}

// Helper function to generate timestamp for filenames
function getTimestampString() {
    const now = new Date();
    const year = now.getFullYear();
    const month = String(now.getMonth() + 1).padStart(2, '0');
    const day = String(now.getDate()).padStart(2, '0');
    const hour = String(now.getHours()).padStart(2, '0');
    const minute = String(now.getMinutes()).padStart(2, '0');
    const second = String(now.getSeconds()).padStart(2, '0');
    return `${year}${month}${day}_${hour}${minute}${second}`;
}

// Export functions
function exportXTImage() {
    if (!xtVisualization || !xtResults) {
        showXTMessage('No simulation results to export.', 'warning');
        return;
    }
    const timestamp = getTimestampString();
    xtVisualization.exportImage(`xtcalc_${timestamp}.png`);
}

function exportXTDataCSV() {
    if (!xtResults) {
        showXTMessage('No simulation results to export.', 'warning');
        return;
    }
    
    // CSV header with all fields
    var csv = 'Time (ms),Position (m),Density (kg/m³),Velocity (m/s),Pressure (kPa),Temperature (K),Gamma,MW (g/mol),Gas ID\n';
    
    for (var j = 0; j < xtResults.xtData.length; j++) {
        var snapshot = xtResults.xtData[j];
        var t = snapshot.t * 1000;
        
        for (var i = 0; i < snapshot.p.length; i++) {
            csv += t.toFixed(4) + ',';                              // Time (ms)
            csv += xtResults.x[i].toFixed(4) + ',';                 // Position (m)
            csv += snapshot.rho[i].toFixed(4) + ',';                // Density (kg/m³)
            csv += snapshot.u[i].toFixed(4) + ',';                  // Velocity (m/s)
            csv += (snapshot.p[i] / 1000).toFixed(4) + ',';         // Pressure (kPa)
            csv += snapshot.T[i].toFixed(4) + ',';                  // Temperature (K)
            csv += snapshot.gamma[i].toFixed(6) + ',';              // Gamma
            csv += snapshot.mw[i].toFixed(4) + ',';                 // Molecular weight (g/mol)
            csv += snapshot.gasId[i] + '\n';                        // Gas ID
        }
    }
    
    const timestamp = getTimestampString();
    downloadXTCSV(csv, `xtcalc_${timestamp}_xtdata.csv`);
}

function exportXTTracerCSV() {
    if (!xtResults) {
        showXTMessage('No simulation results to export.', 'warning');
        return;
    }
    
    var csv = 'Interface,Time (s),Position (m),Velocity (m/s)\n';
    for (var i = 0; i < xtResults.tracers.length; i++) {
        var tracer = xtResults.tracers[i];
        for (var j = 0; j < tracer.trajectory.length; j++) {
            var point = tracer.trajectory[j];
            csv += (i + 1) + ',' + point.t + ',' + point.x + ',' + point.v + '\n';
        }
    }
    
    const timestamp = getTimestampString();
    downloadXTCSV(csv, `xtcalc_${timestamp}_tracers.csv`);
}

function exportXTJSON() {
    if (!xtResults) {
        showXTMessage('No simulation results to export.', 'warning');
        return;
    }
    
    var data = {
        configuration: {
            slabs: xtSlabs,
            gridPoints: parseInt(document.getElementById('xtGridPoints').value),
            simulationTime: parseFloat(document.getElementById('xtSimTime').value),
            cfl: parseFloat(document.getElementById('xtCflNumber').value)
        },
        results: xtResults
    };
    
    var json = JSON.stringify(data, null, 2);
    var blob = new Blob([json], { type: 'application/json' });
    var url = URL.createObjectURL(blob);
    var link = document.createElement('a');
    link.href = url;
    const timestamp = getTimestampString();
    link.download = `xtcalc_${timestamp}.json`;
    link.click();
    URL.revokeObjectURL(url);
}

async function exportXTDataHDF5() {
    if (!xtResults) {
        showXTMessage('No simulation results to export.', 'warning');
        return;
    }
    
    try {
        // Wait for h5wasm to be ready
        if (!window.h5wasm) {
            await window.h5wasmReady;
        }
        
        const { FS } = await window.h5wasm.ready;
        
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
        const f = new window.h5wasm.File(filename, 'w');
        
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
                    const gasId = snapshot.gasId[j];
                    
                    // Create key with tolerance
                    const key = `${gamma.toFixed(6)}_${mw.toFixed(4)}`;
                    
                    if (!speciesMap.has(key)) {
                        speciesMap.set(key, {
                            gamma: gamma,
                            mw: mw,
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
        f.get('metadata').create_dataset({name: 'cfl', data: [parseFloat(document.getElementById('xtCflNumber').value)], dtype: '<f'});
        f.get('metadata').create_dataset({name: 'final_time', data: [xtResults.t], dtype: '<d'});
        f.get('metadata').create_dataset({name: 'integrator', data: [document.getElementById('xtIntegrator').value]});
        f.get('metadata').create_dataset({name: 'interface_method', data: [document.getElementById('xtInterfaceMethod').value]});
        f.get('metadata').create_dataset({name: 'simulation_datetime', data: [now.toISOString()]});
        f.get('metadata').create_dataset({name: 'wall_time', data: [xtResults.wallTime || 0], dtype: '<d'});
        
        // Write species information
        const gammaArray = speciesArray.map(s => s.gamma);
        const mwArray = speciesArray.map(s => s.mw);
        const namesArray = speciesArray.map(s => s.names.join(' / '));
        const gasIdsArray = speciesArray.map(s => s.gasIds[0]);
        
        f.get('species').create_dataset({name: 'gamma', data: gammaArray, dtype: '<d'});
        f.get('species').create_dataset({name: 'molecular_weight', data: mwArray, dtype: '<d'});
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

function downloadXTCSV(csv, filename) {
    var blob = new Blob([csv], { type: 'text/csv' });
    var url = URL.createObjectURL(blob);
    var link = document.createElement('a');
    link.href = url;
    link.download = filename;
    link.click();
    URL.revokeObjectURL(url);
}

function showXTMessage(message, type) {
    var colors = {
        'info': '#0066cc',
        'success': '#28a745',
        'warning': '#ffc107',
        'error': '#dc3545'
    };

    var container = document.getElementById('xtMessageContainer');
    container.innerHTML = '<div style="padding:12px;margin:10px 0;border-radius:3px;background-color:' + colors[type] + '20;border-left:4px solid ' + colors[type] + ';color:' + colors[type] + ';">' + message + '</div>';

    if (type === 'success') {
        setTimeout(function() {
            container.innerHTML = '';
        }, 5000);
    }
}

// ==========================================
// Shocked-Interface Calculation Functions
// ==========================================

let heavyGasMixtureData = [];

function toggleInterfaceSection() {
    const enabled = document.getElementById('enableInterfaceCalc').checked;
    const section = document.getElementById('interfaceSection');
    section.style.display = enabled ? 'block' : 'none';

    if (enabled) {
        // Delay to ensure gasDatabase is loaded
        setTimeout(populateHeavyGasSelect, 0);
    }
}

function populateHeavyGasSelect() {
    const select = document.getElementById('interfaceHeavyGas');

    // Clear existing options except first
    while (select.options.length > 1) {
        select.remove(1);
    }

    // Add standard gases (gasDatabase is an object, not array)
    for (const id in gasDatabase) {
        if (gasDatabase.hasOwnProperty(id)) {
            const gas = gasDatabase[id];
            const option = document.createElement('option');
            option.value = id;
            option.textContent = gas.name;
            select.appendChild(option);
        }
    }

    // Add custom mixture option
    const customOption = document.createElement('option');
    customOption.value = 'custom';
    customOption.textContent = 'Custom Mixture';
    select.appendChild(customOption);
}

function updateInterfaceHeavyGas() {
    const gasId = document.getElementById('interfaceHeavyGas').value;
    const customDiv = document.getElementById('heavyGasCustomMixture');

    if (gasId === 'custom') {
        customDiv.style.display = 'block';
        if (heavyGasMixtureData.length === 0) {
            addHeavyGasMixtureRow();
        }
    } else {
        customDiv.style.display = 'none';
    }
}

function addHeavyGasMixtureRow() {
    const id = 'heavyMix_' + Date.now();
    heavyGasMixtureData.push({ id, gasId: '', fraction: 0 });
    renderHeavyGasMixture();
}

function removeHeavyGasMixtureRow(id) {
    heavyGasMixtureData = heavyGasMixtureData.filter(row => row.id !== id);
    renderHeavyGasMixture();
}

function renderHeavyGasMixture() {
    const container = document.getElementById('heavyGasMixtureList');
    container.innerHTML = '';

    heavyGasMixtureData.forEach((row, index) => {
        const div = document.createElement('div');
        div.style.cssText = 'display: flex; align-items: center; gap: 10px; margin-bottom: 8px;';

        const select = document.createElement('select');
        select.style.width = '150px';
        select.onchange = (e) => { row.gasId = e.target.value; };

        const defaultOpt = document.createElement('option');
        defaultOpt.value = '';
        defaultOpt.textContent = '-- Select Gas --';
        select.appendChild(defaultOpt);

        for (const id in gasDatabase) {
            if (gasDatabase.hasOwnProperty(id)) {
                const gas = gasDatabase[id];
                const opt = document.createElement('option');
                opt.value = id;
                opt.textContent = gas.name;
                opt.selected = id === row.gasId;
                select.appendChild(opt);
            }
        }

        const input = document.createElement('input');
        input.type = 'number';
        input.value = row.fraction;
        input.min = '0';
        input.max = '1';
        input.step = '0.01';
        input.style.width = '80px';
        input.onchange = (e) => { row.fraction = parseFloat(e.target.value); };

        const label = document.createElement('span');
        label.textContent = 'mole fraction';
        label.style.fontSize = '12px';

        const removeBtn = document.createElement('button');
        removeBtn.type = 'button';
        removeBtn.textContent = 'Remove';
        removeBtn.style.cssText = 'font-size: 11px; padding: 4px 8px; background-color: #dc3545; border-color: #dc3545; color: white;';
        removeBtn.onclick = () => removeHeavyGasMixtureRow(row.id);

        div.appendChild(select);
        div.appendChild(input);
        div.appendChild(label);
        div.appendChild(removeBtn);
        container.appendChild(div);
    });
}

function calculateHeavyGasMixture() {
    let totalFraction = 0;
    let avgGamma = 0;
    let avgMW = 0;

    for (const row of heavyGasMixtureData) {
        if (!row.gasId || row.fraction <= 0) {
            alert('Error: All heavy gas mixture components must have a gas selected and positive mole fraction');
            return null;
        }

        const gas = getGasById(row.gasId);
        if (!gas) {
            alert('Error: Invalid gas in heavy gas mixture');
            return null;
        }

        totalFraction += row.fraction;
        avgMW += row.fraction * gas.mw;
        avgGamma += row.fraction * gas.gamma;
    }

    if (Math.abs(totalFraction - 1.0) > 0.01) {
        alert('Error: Heavy gas mole fractions must sum to 1.0 (current sum: ' + totalFraction.toFixed(3) + ')');
        return null;
    }

    return {
        gamma: avgGamma / totalFraction,
        mw: avgMW / totalFraction
    };
}

function calculateInterface() {
    try {
        // Get driven (light) gas properties from main shock calculation
        const drivenGasId = document.Shock.drivengas.value;
        let lightGamma, lightMW;

        if (drivenGasId === 'custom') {
            const drivenProps = calculateMixtureProperties('driven');
            if (!drivenProps) return;
            lightGamma = drivenProps.gamma;
            lightMW = drivenProps.mw;
        } else {
            const drivenGas = getGasById(drivenGasId);
            if (!drivenGas) {
                alert('Error: Please run shock calculation first');
                return;
            }
            lightGamma = drivenGas.gamma;
            lightMW = drivenGas.mw;
        }

        // Get shocked state properties from main calculation
        const p2 = parseFloat(document.Shock.pshocked.value) * 1e6; // Convert MPa to Pa
        const t2 = parseFloat(document.Shock.tshocked.value);
        const u2 = parseFloat(document.Shock.ushocked.value); // Particle velocity in shocked gas

        if (isNaN(p2) || isNaN(t2) || isNaN(u2)) {
            alert('Error: Please run shock calculation first');
            return;
        }

        // Get heavy gas properties
        const heavyGasId = document.getElementById('interfaceHeavyGas').value;
        let heavyGamma, heavyMW;

        if (!heavyGasId) {
            alert('Error: Please select heavy gas');
            return;
        }

        if (heavyGasId === 'custom') {
            const heavyProps = calculateHeavyGasMixture();
            if (!heavyProps) return;
            heavyGamma = heavyProps.gamma;
            heavyMW = heavyProps.mw;
        } else {
            const heavyGas = getGasById(heavyGasId);
            if (!heavyGas) {
                alert('Error: Invalid heavy gas selection');
                return;
            }
            heavyGamma = heavyGas.gamma;
            heavyMW = heavyGas.mw;
        }

        // Get initial conditions for heavy gas
        const p1 = parseFloat(document.Shock.pinit.value);
        const tHeavy = parseFloat(document.getElementById('interfaceHeavyTemp').value);

        if (isNaN(tHeavy) || tHeavy <= 0) {
            alert('Error: Heavy gas temperature must be positive');
            return;
        }

        // Solve Riemann problem at interface
        // Left state: shocked light gas (moving at u2)
        // Right state: ambient heavy gas (at rest)
        const shockedLightState = {
            p: p2,
            T: t2,
            u: u2,  // Shocked gas is already moving
            gamma: lightGamma,
            mw: lightMW
        };

        const ambientHeavyState = {
            p: p1,
            T: tHeavy,
            u: 0,  // Heavy gas initially at rest
            gamma: heavyGamma,
            mw: heavyMW
        };

        const result = solveShockedInterface(shockedLightState, ambientHeavyState);

        // Display results
        document.getElementById('interfaceVelocity').value = result.interfaceVelocity.toPrecision(6);
        document.getElementById('interfacePressure').value = (result.interfacePressure / 1e6).toPrecision(6);
        document.getElementById('atwoodNumber').value = result.atwoodNumber.toPrecision(6);
        document.getElementById('waveConfig').value =
            `Left: ${result.waveConfig.left}, Right: ${result.waveConfig.right}`;

        // Light gas (left) post-shock
        document.getElementById('lightGasPressure').value = (result.lightGasPost.p / 1e6).toPrecision(6);
        document.getElementById('lightGasTemp').value = result.lightGasPost.T.toPrecision(6);
        document.getElementById('lightGasDensity').value = result.lightGasPost.rho.toPrecision(6);
        document.getElementById('lightGasSound').value = result.lightGasPost.a.toPrecision(6);

        // Heavy gas (right) post-shock
        document.getElementById('heavyGasPressure').value = (result.heavyGasPost.p / 1e6).toPrecision(6);
        document.getElementById('heavyGasTemp').value = result.heavyGasPost.T.toPrecision(6);
        document.getElementById('heavyGasDensity').value = result.heavyGasPost.rho.toPrecision(6);
        document.getElementById('heavyGasSound').value = result.heavyGasPost.a.toPrecision(6);

        // Calculate RMI growth predictions
        const lambda = parseFloat(document.getElementById('perturbWavelength').value) / 1000; // mm to m
        const a0 = parseFloat(document.getElementById('perturbAmplitude').value) / 1000; // mm to m
        const downstreamDist = parseFloat(document.getElementById('rmiDownstreamDist').value) / 1000; // mm to m

        const deltaU = calculateVelocityJump(result.interfaceVelocity, 0);

        // Calculate time to reach downstream distance
        const tMax = downstreamDist / result.interfaceVelocity;

        const rmiParams = calculateRMIParameters({
            eta0: a0,
            lambda: lambda,
            A: result.atwoodNumber,
            deltaU: deltaU
        });

        // Calculate amplitude at downstream location using each model
        const k = 2 * Math.PI / lambda;
        const amplitudeLinear = rmiLinearAmplitude(a0, result.atwoodNumber, k, deltaU, tMax);
        // const amplitudeZS = rmiZhangSohnAmplitude(a0, result.atwoodNumber, k, deltaU, tMax);
        const amplitudeDR = rmiDimonteAmplitude(a0, result.atwoodNumber, k, deltaU, tMax);

        // Display RMI parameters
        document.getElementById('rmiWavenumber').textContent = rmiParams.wavenumber.toPrecision(4);
        document.getElementById('rmiKa0').textContent = rmiParams.ka0.toPrecision(4);
        document.getElementById('rmiGrowthRate').textContent = rmiParams.initialGrowthRate.toPrecision(4);
        document.getElementById('rmiDoublingTime').textContent = (rmiParams.doublingTime * 1000).toPrecision(4);
        document.getElementById('rmiTimeAtLocation').textContent = (tMax * 1000).toPrecision(4);
        document.getElementById('rmiAmplitudeAtLocation').textContent =
            `Linear: ${(amplitudeLinear * 1000).toPrecision(3)} | DR: ${(amplitudeDR * 1000).toPrecision(3)}`;

        // Generate and plot growth curves
        const growth = calculateGrowthTimeSeries({
            eta0: a0,
            lambda: lambda,
            A: result.atwoodNumber,
            deltaU: deltaU,
            tMax: tMax,
            numPoints: 200
        });

        plotRMIGrowth(growth, downstreamDist, result.interfaceVelocity);

        // Show results section
        document.getElementById('interfaceResults').style.display = 'block';

    } catch (error) {
        console.error('Error in interface calculation:', error);
        alert('An error occurred during interface calculation: ' + error.message);
    }
}

function plotRMIGrowth(growth, downstreamDist, interfaceVelocity) {
    const canvas = document.getElementById('rmiGrowthCanvas');
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;

    // Clear canvas
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, width, height);

    // Convert times to ms and amplitudes to mm for plotting
    const timesMs = growth.times.map(t => t * 1000);
    const linearMm = growth.linear.map(a => a * 1000);
    // const zhangSohnMm = growth.zhangSohn.map(a => a * 1000);
    const dimonteMm = growth.dimonte.map(a => a * 1000);

    // Time at downstream location
    const tDownstream = (downstreamDist / interfaceVelocity) * 1000; // ms

    // Find data ranges
    const tMin = 0;
    const tMax = Math.max(...timesMs);
    const aMin = Math.min(0, ...linearMm, ...dimonteMm);
    const aMax = Math.max(...linearMm, ...dimonteMm);

    // Add 10% padding
    const aRange = aMax - aMin;
    const aMinPlot = aMin - 0.1 * aRange;
    const aMaxPlot = aMax + 0.1 * aRange;

    // Margins
    const margin = { left: 80, right: 40, top: 40, bottom: 60 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    // Scaling functions
    const scaleX = (t) => margin.left + (t - tMin) / (tMax - tMin) * plotWidth;
    const scaleY = (a) => margin.top + plotHeight - (a - aMinPlot) / (aMaxPlot - aMinPlot) * plotHeight;

    // Draw axes
    ctx.strokeStyle = '#000000';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(margin.left, margin.top);
    ctx.lineTo(margin.left, margin.top + plotHeight);
    ctx.lineTo(margin.left + plotWidth, margin.top + plotHeight);
    ctx.stroke();

    // Draw grid
    ctx.strokeStyle = '#e0e0e0';
    ctx.lineWidth = 1;
    const numGridLines = 5;

    for (let i = 0; i <= numGridLines; i++) {
        const t = tMin + (tMax - tMin) * i / numGridLines;
        const x = scaleX(t);
        ctx.beginPath();
        ctx.moveTo(x, margin.top);
        ctx.lineTo(x, margin.top + plotHeight);
        ctx.stroke();

        const a = aMinPlot + (aMaxPlot - aMinPlot) * i / numGridLines;
        const y = scaleY(a);
        ctx.beginPath();
        ctx.moveTo(margin.left, y);
        ctx.lineTo(margin.left + plotWidth, y);
        ctx.stroke();
    }

    // Plot curves
    function plotCurve(data, color, lineWidth = 2) {
        ctx.strokeStyle = color;
        ctx.lineWidth = lineWidth;
        ctx.beginPath();
        for (let i = 0; i < timesMs.length; i++) {
            const x = scaleX(timesMs[i]);
            const y = scaleY(data[i]);
            if (i === 0) {
                ctx.moveTo(x, y);
            } else {
                ctx.lineTo(x, y);
            }
        }
        ctx.stroke();
    }

    plotCurve(linearMm, '#0066cc', 2);
    // plotCurve(zhangSohnMm, '#dc3545', 2);
    plotCurve(dimonteMm, '#28a745', 2);

    // Mark downstream location with vertical line
    if (tDownstream >= tMin && tDownstream <= tMax) {
        const xDown = scaleX(tDownstream);
        ctx.strokeStyle = '#ff9800';
        ctx.lineWidth = 2;
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(xDown, margin.top);
        ctx.lineTo(xDown, margin.top + plotHeight);
        ctx.stroke();
        ctx.setLineDash([]);

        // Label
        ctx.fillStyle = '#ff9800';
        ctx.font = 'bold 11px Arial';
        ctx.textAlign = 'center';
        ctx.fillText('Downstream', xDown, margin.top - 5);
    }

    // Axis labels
    ctx.fillStyle = '#000000';
    ctx.font = '14px Arial';
    ctx.textAlign = 'center';
    ctx.fillText('Time (ms)', margin.left + plotWidth / 2, height - 10);

    ctx.save();
    ctx.translate(15, margin.top + plotHeight / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('Amplitude (mm)', 0, 0);
    ctx.restore();

    // Title
    ctx.font = 'bold 16px Arial';
    ctx.fillText('RMI Perturbation Growth', width / 2, 25);

    // Tick labels
    ctx.font = '12px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (let i = 0; i <= numGridLines; i++) {
        const t = tMin + (tMax - tMin) * i / numGridLines;
        const x = scaleX(t);
        ctx.fillText(t.toFixed(1), x, margin.top + plotHeight + 5);
    }

    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let i = 0; i <= numGridLines; i++) {
        const a = aMinPlot + (aMaxPlot - aMinPlot) * i / numGridLines;
        const y = scaleY(a);
        ctx.fillText(a.toFixed(2), margin.left - 10, y);
    }
}

