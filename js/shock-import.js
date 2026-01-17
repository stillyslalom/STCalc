/**
 * Import shock tube configuration from the existing calculator
 * Reads from localStorage and provides gas database
 */

// Gas database with properties (matches the main calculator)
const GasDatabase = {
    "1": { name: "Air", gamma: 1.4020, mw: 28.9660 },
    "2": { name: "Argon", gamma: 1.6700, mw: 39.948000 },
    "3": { name: "Carbon Dioxide", gamma: 1.29700, mw: 44.01000 },
    "4": { name: "Helium", gamma: 1.66700, mw: 4.002602000 },
    "5": { name: "Neon", gamma: 1.66700, mw: 20.179700 },
    "6": { name: "Nitrogen", gamma: 1.40100, mw: 28.0134000 },
    "7": { name: "Xenon", gamma: 1.66700, mw: 131.29000 },
    "8": { name: "Sulfur Hexafluoride", gamma: 1.09200, mw: 146.06000 },
    "9": { name: "Helium-Neon", gamma: 1.668, mw: 21.98 },
    "10": { name: "Acetone Vapor", gamma: 1.14000, mw: 58.08000 }
};

/**
 * Load custom gases from localStorage
 */
function loadCustomGases() {
    try {
        const stored = localStorage.getItem('wistl_custom_gases');
        if (stored) {
            const gases = JSON.parse(stored);
            const customGases = {};
            for (let gas of gases) {
                customGases[gas.id] = {
                    name: gas.name,
                    gamma: gas.gamma,
                    mw: gas.mw
                };
            }
            return customGases;
        }
    } catch (e) {
        console.error('Error loading custom gases:', e);
    }
    return {};
}

/**
 * Get all available gases (predefined + custom)
 */
function getAllGases() {
    const allGases = { ...GasDatabase };
    const customGases = loadCustomGases();
    return { ...allGases, ...customGases };
}

/**
 * Calculate mixture properties from components
 */
function calculateMixtureProperties(components) {
    if (!components || components.length === 0) {
        return null;
    }
    
    const allGases = getAllGases();
    
    // Validate total fraction
    const totalFraction = components.reduce((sum, c) => sum + c.fraction, 0);
    if (Math.abs(totalFraction - 1.0) > 0.001) {
        console.error('Mixture fractions do not sum to 1.0');
        return null;
    }
    
    // Calculate mixture molecular weight
    let mw_mix = 0;
    for (let component of components) {
        const gas = allGases[component.gasId];
        if (!gas) {
            console.error(`Gas ${component.gasId} not found`);
            return null;
        }
        mw_mix += component.fraction * gas.mw;
    }
    
    // Calculate mixture gamma using specific heats
    let cv_mix = 0;
    let cp_mix = 0;
    
    for (let component of components) {
        const gas = allGases[component.gasId];
        const x = component.fraction;
        const gamma = gas.gamma;
        
        const cv = 1.0 / (gamma - 1.0);
        const cp = gamma * cv;
        
        cv_mix += x * cv;
        cp_mix += x * cp;
    }
    
    const gamma_mix = cp_mix / cv_mix;
    
    return {
        gamma: gamma_mix,
        mw: mw_mix
    };
}

/**
 * Import calculator state from localStorage
 */
function importFromCalculator() {
    try {
        const stored = localStorage.getItem('wistl_calculator_state');
        if (!stored) {
            return null;
        }
        
        const state = JSON.parse(stored);
        const allGases = getAllGases();
        
        // Process driven gas
        let drivenGas = null;
        if (state.drivengas === 'custom' && state.drivenMixture) {
            const props = calculateMixtureProperties(state.drivenMixture);
            if (props) {
                drivenGas = {
                    gasId: 'custom_driven',
                    name: 'Custom Mixture (Driven)',
                    gamma: props.gamma,
                    mw: props.mw
                };
            }
        } else if (state.drivengas && allGases[state.drivengas]) {
            const gas = allGases[state.drivengas];
            drivenGas = {
                gasId: state.drivengas,
                name: gas.name,
                gamma: gas.gamma,
                mw: gas.mw
            };
        }
        
        // Process driver gas
        let driverGas = null;
        if (state.drivergas === 'custom' && state.driverMixture) {
            const props = calculateMixtureProperties(state.driverMixture);
            if (props) {
                driverGas = {
                    gasId: 'custom_driver',
                    name: 'Custom Mixture (Driver)',
                    gamma: props.gamma,
                    mw: props.mw
                };
            }
        } else if (state.drivergas && allGases[state.drivergas]) {
            const gas = allGases[state.drivergas];
            driverGas = {
                gasId: state.drivergas,
                name: gas.name,
                gamma: gas.gamma,
                mw: gas.mw
            };
        }
        
        if (!drivenGas || !driverGas) {
            console.error('Failed to parse gas properties');
            return null;
        }
        
        return {
            driven: {
                gas: drivenGas,
                pressure: parseFloat(state.pinit) || 101325,
                temperature: parseFloat(state.tinit) || 300
            },
            driver: {
                gas: driverGas,
                pressure: parseFloat(state.pfour) * 1e6 || 300000,  // pfour is in MPa, convert to Pa
                temperature: parseFloat(state.tid) || 300
            },
            mach: parseFloat(state.mach) || 1.5
        };
    } catch (e) {
        console.error('Error importing calculator state:', e);
        return null;
    }
}

/**
 * Populate gas dropdown with all available gases
 */
function populateGasDropdown(selectElement) {
    const allGases = getAllGases();
    
    // Clear existing options
    selectElement.innerHTML = '';
    
    // Add empty option
    const emptyOption = document.createElement('option');
    emptyOption.value = '';
    emptyOption.textContent = 'Select Gas';
    selectElement.appendChild(emptyOption);
    
    // Sort gases alphabetically
    const gasEntries = Object.entries(allGases).sort((a, b) => 
        a[1].name.localeCompare(b[1].name)
    );
    
    // Add gas options
    for (let [id, gas] of gasEntries) {
        const option = document.createElement('option');
        option.value = id;
        option.textContent = id.startsWith('custom_') ? `★ ${gas.name}` : gas.name;
        selectElement.appendChild(option);
    }
}

/**
 * Get gas properties by ID
 */
function getGasById(gasId) {
    const allGases = getAllGases();
    return allGases[gasId] || null;
}
