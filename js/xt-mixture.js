/**
 * X-T Diagram Mixture Support
 * Handles custom gas mixtures for X-T diagram slabs
 */

// Mixture counter for unique IDs
var xtMixtureComponentCounter = 0;

/**
 * Toggle between single gas and custom mixture for an X-T slab
 */
function toggleXTSlabMixture(slabId) {
    var slab = xtSlabs.find(s => s.id === slabId);
    if (!slab) return;
    
    var select = document.getElementById('xt_slab_gas_' + slabId);
    if (!select) return;
    
    var mixtureDiv = document.getElementById('xt_slab_mixture_' + slabId);
    
    if (select.value === 'custom_mixture') {
        // Show mixture inputs
        if (!mixtureDiv) {
            mixtureDiv = document.createElement('div');
            mixtureDiv.id = 'xt_slab_mixture_' + slabId;
            mixtureDiv.className = 'custom-gas-section';
            mixtureDiv.style.marginTop = '10px';
            mixtureDiv.style.gridColumn = '1 / -1'; // Span all columns
            
            var slabItem = document.getElementById('xt_slab_item_' + slabId);
            slabItem.appendChild(mixtureDiv);
        }
        
        // Initialize mixture if needed
        if (!xtSlabMixtures[slabId]) {
            xtSlabMixtures[slabId] = [];
            addXTMixtureComponent(slabId);
        }
        
        renderXTMixtureInputs(slabId);
        mixtureDiv.style.display = 'block';
    } else {
        // Hide mixture inputs
        if (mixtureDiv) {
            mixtureDiv.style.display = 'none';
        }
    }
}

/**
 * Add a new component to an X-T slab mixture
 */
function addXTMixtureComponent(slabId) {
    if (!xtSlabMixtures[slabId]) {
        xtSlabMixtures[slabId] = [];
    }
    
    var compId = xtMixtureComponentCounter++;
    xtSlabMixtures[slabId].push({
        id: compId,
        gasId: '',
        fraction: 0
    });
    
    renderXTMixtureInputs(slabId);
}

/**
 * Remove a component from an X-T slab mixture
 */
function removeXTMixtureComponent(slabId, compId) {
    if (!xtSlabMixtures[slabId]) return;
    
    xtSlabMixtures[slabId] = xtSlabMixtures[slabId].filter(c => c.id !== compId);
    renderXTMixtureInputs(slabId);
}

/**
 * Update a mixture component property
 */
function updateXTMixtureComponent(slabId, compId, property, value) {
    if (!xtSlabMixtures[slabId]) return;
    
    var comp = xtSlabMixtures[slabId].find(c => c.id === compId);
    if (comp) {
        comp[property] = value;
        updateXTMixtureTotal(slabId);
    }
}

/**
 * Handle gas selection for X-T mixture component
 */
function handleXTMixtureGasSelection(slabId, compId, value) {
    if (value === 'define_new') {
        showXTDefineGasForm(slabId, compId);
    } else {
        // Hide definition form if showing
        var defineDiv = document.getElementById('xt_define_' + slabId + '_' + compId);
        if (defineDiv) {
            defineDiv.style.display = 'none';
        }
        updateXTMixtureComponent(slabId, compId, 'gasId', value);
    }
}

/**
 * Show custom gas definition form for X-T mixture
 */
function showXTDefineGasForm(slabId, compId) {
    var defineDiv = document.getElementById('xt_define_' + slabId + '_' + compId);
    if (!defineDiv) return;
    
    defineDiv.innerHTML = '<div style="background-color:#fff3cd;border:1px solid #ffc107;border-radius:3px;padding:10px;margin-top:5px;">' +
        '<b>Define New Gas:</b><br>' +
        'Name: <input type="text" id="xt_newgas_name_' + slabId + '_' + compId + '" size="15" placeholder="e.g., Methane" style="margin:5px 5px 5px 0;"> ' +
        'Gamma (γ): <input type="text" id="xt_newgas_gamma_' + slabId + '_' + compId + '" size="8" placeholder="e.g., 1.31" style="margin:5px 5px 5px 0;"> ' +
        'MW (g/mol): <input type="text" id="xt_newgas_mw_' + slabId + '_' + compId + '" size="8" placeholder="e.g., 16.04" style="margin:5px 5px 5px 0;"> ' +
        '<button type="button" onclick="saveXTCustomGas(' + slabId + ', ' + compId + ')" style="margin:5px 5px 5px 0;">Save Gas</button> ' +
        '<button type="button" onclick="cancelXTDefineGas(' + slabId + ', ' + compId + ')" style="margin:5px 0;">Cancel</button>' +
        '</div>';
    
    defineDiv.style.display = 'block';
}

/**
 * Save a new custom gas from X-T mixture interface
 */
function saveXTCustomGas(slabId, compId) {
    var nameInput = document.getElementById('xt_newgas_name_' + slabId + '_' + compId);
    var gammaInput = document.getElementById('xt_newgas_gamma_' + slabId + '_' + compId);
    var mwInput = document.getElementById('xt_newgas_mw_' + slabId + '_' + compId);
    
    var name = nameInput.value.trim();
    var gamma = parseFloat(gammaInput.value);
    var mw = parseFloat(mwInput.value);
    
    // Validate inputs
    if (!name) {
        alert('Please enter a gas name.');
        nameInput.focus();
        return;
    }
    
    if (isNaN(gamma) || gamma <= 1.0) {
        alert('Gamma must be a number greater than 1.0 (thermodynamically required for ideal gases).');
        gammaInput.focus();
        return;
    }
    
    if (isNaN(mw) || mw <= 0) {
        alert('Molecular weight must be a positive number.');
        mwInput.focus();
        return;
    }
    
    // Check for duplicate name
    var allGases = getAllGases();
    for (var id in allGases) {
        if (allGases[id].name.toLowerCase() === name.toLowerCase()) {
            alert('A gas with this name already exists. Please use a different name.');
            nameInput.focus();
            return;
        }
    }
    
    // Create new custom gas
    var gasId = 'custom_' + customGasCounter++;
    customGases[gasId] = {
        id: gasId,
        name: name,
        gamma: gamma,
        mw: mw
    };
    
    // Save to localStorage
    saveCustomGases();
    
    // Update the component's gas selection
    updateXTMixtureComponent(slabId, compId, 'gasId', gasId);
    
    // Hide the definition form
    var defineDiv = document.getElementById('xt_define_' + slabId + '_' + compId);
    if (defineDiv) {
        defineDiv.style.display = 'none';
    }
    
    // Refresh all X-T mixture dropdowns to include the new gas
    refreshAllXTMixtureDropdowns();
    
    // Update custom gases list in main calculator if present
    if (typeof updateCustomGasList === 'function') {
        updateCustomGasList();
    }
    
    // Re-render to show the selected gas
    renderXTMixtureInputs(slabId);
    
    alert('Custom gas "' + name + '" saved successfully!');
}

/**
 * Cancel custom gas definition in X-T mixture
 */
function cancelXTDefineGas(slabId, compId) {
    var select = document.getElementById('xt_mix_gas_' + slabId + '_' + compId);
    var defineDiv = document.getElementById('xt_define_' + slabId + '_' + compId);
    
    if (select) {
        select.value = '';
    }
    if (defineDiv) {
        defineDiv.style.display = 'none';
    }
    
    updateXTMixtureComponent(slabId, compId, 'gasId', '');
}

/**
 * Refresh all X-T mixture dropdowns to include newly added custom gases
 */
function refreshAllXTMixtureDropdowns() {
    // Iterate through all slabs with custom mixtures
    for (var slabId in xtSlabMixtures) {
        var components = xtSlabMixtures[slabId];
        for (var i = 0; i < components.length; i++) {
            var comp = components[i];
            var select = document.getElementById('xt_mix_gas_' + slabId + '_' + comp.id);
            if (select) {
                var currentValue = select.value;
                
                // Rebuild options
                var allGases = getAllGases();
                var html = '<option value="">Select Gas</option>';
                html += '<option value="define_new">--- Define New Gas... ---</option>';
                html += '<option disabled>──────────</option>';
                
                var gasEntries = Object.entries(allGases).sort(function(a, b) {
                    return a[1].name.localeCompare(b[1].name);
                });
                
                for (var j = 0; j < gasEntries.length; j++) {
                    var gasId = gasEntries[j][0];
                    var gas = gasEntries[j][1];
                    var prefix = gasId.startsWith('custom_') ? '★ ' : '';
                    html += '<option value="' + gasId + '">' + prefix + gas.name + '</option>';
                }
                
                select.innerHTML = html;
                select.value = currentValue;
            }
        }
    }
}

/**
 * Render mixture input controls for a slab
 */
function renderXTMixtureInputs(slabId) {
    var mixtureDiv = document.getElementById('xt_slab_mixture_' + slabId);
    if (!mixtureDiv) return;
    
    var components = xtSlabMixtures[slabId] || [];
    
    var html = '<div style="font-weight: 600; margin-bottom: 8px; font-size: 13px;">Custom Mixture (Mole Fractions):</div>';
    
    // Render each component
    for (var i = 0; i < components.length; i++) {
        var comp = components[i];
        html += '<div class="mixture-component" style="display: flex; flex-wrap: nowrap; align-items: center; gap: 8px; margin-bottom: 8px;">';
        html += '<select id="xt_mix_gas_' + slabId + '_' + comp.id + '" ';
        html += 'onchange="handleXTMixtureGasSelection(' + slabId + ', ' + comp.id + ', this.value)" ';
        html += 'style="min-width: 150px;">';
        html += '<option value="">Select Gas</option>';
        html += '<option value="define_new">--- Define New Gas... ---</option>';
        html += '<option disabled>──────────</option>';
        
        // Add all available gases
        var allGases = getAllGases();
        var gasEntries = Object.entries(allGases).sort(function(a, b) {
            return a[1].name.localeCompare(b[1].name);
        });
        
        for (var j = 0; j < gasEntries.length; j++) {
            var gasId = gasEntries[j][0];
            var gas = gasEntries[j][1];
            var selected = gasId === comp.gasId ? 'selected' : '';
            var prefix = gasId.startsWith('custom_') ? '★ ' : '';
            html += '<option value="' + gasId + '" ' + selected + '>' + prefix + gas.name + '</option>';
        }
        
        html += '</select>';
        html += '<span style="font-size: 13px; white-space: nowrap;">Mole Fraction:</span>';
        html += '<input type="number" id="xt_mix_frac_' + slabId + '_' + comp.id + '" ';
        html += 'value="' + (comp.fraction || 0) + '" ';
        html += 'onchange="updateXTMixtureComponent(' + slabId + ', ' + comp.id + ', \'fraction\', parseFloat(this.value))" ';
        html += 'step="0.01" min="0" max="1" style="width: 80px;">';
        html += '<button type="button" onclick="removeXTMixtureComponent(' + slabId + ', ' + comp.id + ')" ';
        html += 'style="padding: 4px 10px; font-size: 12px; white-space: nowrap;">Remove</button>';
        
        // Container for custom gas definition form
        html += '<div id="xt_define_' + slabId + '_' + comp.id + '" style="display:none;"></div>';
        
        html += '</div>';
    }
    
    // Add component button
    html += '<button type="button" onclick="addXTMixtureComponent(' + slabId + ')" ';
    html += 'style="margin-top: 8px; padding: 6px 12px; font-size: 13px;">Add Gas Component</button>';
    
    // Total display
    html += '<span id="xt_mix_total_' + slabId + '" style="margin-left: 10px; font-size: 13px;"></span>';
    
    mixtureDiv.innerHTML = html;
    updateXTMixtureTotal(slabId);
}

/**
 * Update and display the total mole fraction for a slab mixture
 */
function updateXTMixtureTotal(slabId) {
    var totalSpan = document.getElementById('xt_mix_total_' + slabId);
    if (!totalSpan) return;
    
    var components = xtSlabMixtures[slabId] || [];
    var total = 0;
    
    for (var i = 0; i < components.length; i++) {
        var frac = parseFloat(components[i].fraction) || 0;
        total += frac;
    }
    
    var displayTotal = Math.round(total * 100 * 100) / 100;
    
    if (Math.abs(total - 1.0) < 0.001) {
        totalSpan.innerHTML = '<span style="color:#28a745;">Total: ' + displayTotal + '% ✓</span>';
    } else if (total > 0) {
        totalSpan.innerHTML = '<span style="color:#dc3545;">Total: ' + displayTotal + '% (must equal 100%)</span>';
    } else {
        totalSpan.innerHTML = '<span style="color:#6c757d;">Total: 0%</span>';
    }
}

/**
 * Calculate effective properties for a custom mixture slab
 */
function calculateXTMixtureProperties(slabId) {
    var components = xtSlabMixtures[slabId];
    if (!components || components.length === 0) {
        return null;
    }
    
    var totalFraction = 0;
    var validComponents = [];
    var allGases = getAllGases();
    
    // Collect valid components
    for (var i = 0; i < components.length; i++) {
        var comp = components[i];
        var fraction = parseFloat(comp.fraction) || 0;
        
        if (comp.gasId && fraction > 0) {
            var gas = allGases[comp.gasId];
            if (gas) {
                validComponents.push({
                    gas: gas,
                    fraction: fraction
                });
                totalFraction += fraction;
            }
        }
    }
    
    // Validate total
    if (Math.abs(totalFraction - 1.0) > 0.001) {
        throw new Error('Mixture mole fractions must sum to 1.0 (100%). Current total: ' + 
                       (totalFraction * 100).toFixed(2) + '%');
    }
    
    if (validComponents.length === 0) {
        throw new Error('No valid gas components specified in mixture.');
    }
    
    // Calculate mixture molecular weight
    var mw_mix = 0;
    for (var i = 0; i < validComponents.length; i++) {
        mw_mix += validComponents[i].fraction * validComponents[i].gas.mw;
    }
    
    // Calculate mixture gamma using specific heats
    var cv_mix = 0;
    var cp_mix = 0;
    
    for (var i = 0; i < validComponents.length; i++) {
        var x = validComponents[i].fraction;
        var gamma = validComponents[i].gas.gamma;
        
        var cv = 1.0 / (gamma - 1.0);
        var cp = gamma * cv;
        
        cv_mix += x * cv;
        cp_mix += x * cp;
    }
    
    var gamma_mix = cp_mix / cv_mix;
    
    return {
        gamma: gamma_mix,
        mw: mw_mix,
        components: validComponents
    };
}

/**
 * Get gas properties for a slab (either single gas or mixture)
 */
function getXTSlabGasProperties(slab) {
    if (slab.gas === 'custom_mixture') {
        // Calculate mixture properties
        var props = calculateXTMixtureProperties(slab.id);
        
        // Generate human-readable name
        var allGases = getAllGases();
        var names = [];
        for (var i = 0; i < props.components.length; i++) {
            var comp = props.components[i];
            var percent = Math.round(comp.fraction * 100);
            names.push(percent + '% ' + comp.gas.name);
        }
        
        return {
            gamma: props.gamma,
            mw: props.mw,
            gasName: names.join(' - ')
        };
    } else {
        // Single gas
        var gas = getAllGases()[slab.gas];
        if (!gas) {
            throw new Error('Invalid gas selection');
        }
        return {
            gamma: gas.gamma,
            mw: gas.mw,
            gasName: gas.name
        };
    }
}

/**
 * Get display name for a slab's gas
 */
function getXTSlabGasDisplayName(slab) {
    if (slab.gas === 'custom_mixture') {
        var components = xtSlabMixtures[slab.id];
        if (!components || components.length === 0) {
            return 'Custom Mixture (not configured)';
        }
        
        var allGases = getAllGases();
        var names = [];
        
        for (var i = 0; i < components.length; i++) {
            var comp = components[i];
            if (comp.gasId && comp.fraction > 0) {
                var gas = allGases[comp.gasId];
                if (gas) {
                    var percent = Math.round(comp.fraction * 100);
                    names.push(percent + '% ' + gas.name);
                }
            }
        }
        
        return 'Custom Mixture: ' + (names.length > 0 ? names.join(' + ') : 'not configured');
    } else {
        var gas = getAllGases()[slab.gas];
        return gas ? gas.name : 'Not selected';
    }
}
