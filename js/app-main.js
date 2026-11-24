


function calcShock() {
    try {
        var p1,p2,p5,p2p1,p5p2;  //pressure and ratio

        var t1,t2,t5,t2t1,t5t2;  //temperature and ratio

        var d1,d2,d5,d2d1,d5d2;  //density and ratio

        var a1,a2,a5;            //speed of sound

        var m,g,gp,gm,mw,ru,r;

        var u2,w;

        var mr,wr;

        var g4,mw4,a4,p4p1,p4,r4,t4,areafactor;

        var d4,p4p1exp,p4p1arg;

        m  = parseFloat(document.Shock.mach.value);

        p1 = parseFloat(document.Shock.pinit.value);

        t1 = parseFloat(document.Shock.tinit.value);

        areafactor = parseFloat(document.Shock.areafact.value);

        // Validate input parameters for physical constraints
        if (isNaN(m) || m <= 0) {
            alert('Error: Mach number must be a positive number greater than 0');
            return;
        }

        if (isNaN(p1) || p1 <= 0) {
            alert('Error: Initial pressure must be a positive number greater than 0');
            return;
        }

        if (isNaN(t1) || t1 <= 0) {
            alert('Error: Initial temperature must be a positive number greater than 0 K');
            return;
        }

        if (isNaN(areafactor) || areafactor <= 0) {
            alert('Error: Area factor must be a positive number greater than 0');
            return;
        }

        // Get driven gas properties
        var drivenGasId = document.Shock.drivengas.value;

        if (drivenGasId === 'custom') {
            // Calculate properties for custom mixture
            var drivenProps = calculateMixtureProperties('driven');
            if (drivenProps === null) {
                return; // Error already shown in calculateMixtureProperties
            }
            g = drivenProps.gamma;
            mw = drivenProps.mw;
        } else {
            // Use gas database lookup
            var drivenGas = getGasById(drivenGasId);
            if (!drivenGas) {
                alert('Error: Invalid driven gas selection');
                return;
            }
            g = drivenGas.gamma;
            mw = drivenGas.mw;
        }

        // Get driver gas properties
        var driverGasId = document.Shock.drivergas.value;

        if (driverGasId === 'custom') {
            // Calculate properties for custom mixture
            var driverProps = calculateMixtureProperties('driver');
            if (driverProps === null) {
                return; // Error already shown in calculateMixtureProperties
            }
            g4 = driverProps.gamma;
            mw4 = driverProps.mw;
        } else {
            // Use gas database lookup
            var driverGas = getGasById(driverGasId);
            if (!driverGas) {
                alert('Error: Invalid driver gas selection');
                return;
            }
            g4 = driverGas.gamma;
            mw4 = driverGas.mw;
        }

      	t4 = parseFloat(document.Shock.tid.value);

        // Validate driver temperature
        if (isNaN(t4) || t4 <= 0) {
            alert('Error: Driver temperature must be a positive number greater than 0 K');
            return;
        }

	//document.write("t4 " + t4 + "  ");    

        gp = g + 1.00;

        gm = g - 1.00;

      	//document.write("gm " + gm + "  ");

        //document.write("g " + g + "  ");

        //document.write("mw " + mw + "  ");

        ru=parseFloat(8314.510);

        //document.write("ru " + ru + "  ");

        r=ru/mw;

        //document.write("r " + r + "  ");

        a1=Math.sqrt(g*r*t1);

        w=a1*m;

        d1=p1/r/t1;

        //document.write("p2 " + t1 + "  ");

        //document.write("gp " + gp + "  ");

        //document.write("gm " + gm + "  ");

        //ees p2p1=1+((2*gamma_1)/gp_1)*(M_init^2-1);

        //

        //  State 2 properties

        //

        p2p1=1+2*g/gp*(Math.pow(m,2)-1);

        p2=p2p1*p1;

        t2t1=p2p1*((gp/gm+p2p1)/(1+gp/gm*p2p1));

        t2=t2t1*t1;

        a2=Math.sqrt(g*r*t2);

        //ees rho2rho1=(1+gp_1/gm_1*p2p1)/(gp_1/gm_1+p2p1)

	d2d1=(1+(gp/gm)*p2p1)/(gp/gm+p2p1);

	d2=d2d1*d1;

	//ees u_2=c_1/gamma_1*(p2p1-1)*((2*gamma_1/gp_1)/(p2p1+gm_1/gp_1))^0.5

	u2=a1/g*(p2p1-1)*(Math.sqrt((2*g/gp)/(p2p1+gm/gp)));

        //

        //  State 4 properties

        //

	//document.write("gm " + gm + "  ");

        r4=ru/mw4;

        a4=Math.sqrt(g4*r4*t4);

        //ees p4p1=p2p1*(1-((gamma_4-1)*(a_1/a_4)*(p2p1-1))/sqrt(2*gamma_1*(2*gamma_1+;

	//         (gamma_1+1)*(p2p1-1))))^(-2*gamma_4/(gamma_4-1)) ;

        p4p1exp=(-2*g4/(g4-1));

        p4p1arg=(1-((g4-1)*(a1/a4)*(p2p1-1))/Math.sqrt(2*g*(2*g+gp*(p2p1-1))));

	p4p1=p2p1*Math.pow(p4p1arg,p4p1exp);

	p4=p4p1*p1;

	d4=p4/r4/t4;

        //document.write("t4 " + t4 + "  ");

        //document.write("a4 " + a4 + "  ");

        //

        //  State 5 properties

        //

        p5p2=((3*g-1)*p2p1-gm)/(gm*p2p1+gp);

        p5=p5p2*p2;

        t5t2=p5p2*((gp/gm+p5p2)/(1+gp/gm*p5p2));

        t5=t5t2*t2;

        d5d2=(1+(gp/gm)*p5p2)/(gp/gm+p5p2);

	d5=d5d2*d2;

        a5=Math.sqrt(g*r*t5);

        wr=(d2*u2)/(d5-d2);

        mr=(wr+u2)/a2;

        document.Shock.pdriv.value = (p1/1e6).toPrecision(6);

        document.Shock.tdriv.value = t1.toPrecision(6);

        document.Shock.ddriv.value = d1.toPrecision(6);

        document.Shock.adriv.value = a1.toPrecision(6);

        document.Shock.pshocked.value = (p2/1e6).toPrecision(6);

        document.Shock.tshocked.value = t2.toPrecision(6);

        document.Shock.dshocked.value = d2.toPrecision(6);

        document.Shock.ashocked.value = a2.toPrecision(6);

        document.Shock.ushocked.value = u2.toPrecision(6);

        document.Shock.pfour.value = (p4/1e6).toPrecision(6);

        document.Shock.pfourpsig.value = ((((p4-101325.0)/101325.0)*14.696)/areafactor).toPrecision(6);

        document.Shock.tfour.value = t4.toPrecision(6);

        document.Shock.dfour.value = d4.toPrecision(6);

        document.Shock.afour.value = a4.toPrecision(6);

        document.Shock.preflect.value = (p5/1e6).toPrecision(6);

        document.Shock.treflect.value = t5.toPrecision(6);

        document.Shock.dreflect.value = d5.toPrecision(6);

        document.Shock.areflect.value = a5.toPrecision(6);

        document.Shock.wreflect.value = wr.toPrecision(6);

        document.Shock.mreflect.value = mr.toPrecision(6);

        document.Shock.wshocked.value = w.toPrecision(6);

        document.Shock.outmw.value = mw.toPrecision(6);

        document.Shock.outgamma.value = g.toPrecision(6);

        // Save calculator state to localStorage
        saveCalculatorState();

    } catch (error) {
        console.error('Error in shock calculation:', error);
        alert('An error occurred during shock tube calculation: ' + error.message +
              '\n\nPlease check your input values and try again.');
    }
}



function calcTime() {
    try {
        var xwave,xtrig,twave,wactual,mactual,xint;

        var tdelayshock,tdelaymembrane,u2;

        xwave = parseFloat(document.Timing.xwave.value);

        twave = parseFloat(document.Timing.twave.value)/1000.00;

        u2 = parseFloat(document.Shock.ushocked.value);

        wactual=xwave/twave;

        mactual=wactual/parseFloat(document.Shock.adriv.value);

        xtrig = parseFloat(document.Timing.xtrig.value);

        xint = parseFloat(document.Timing.xint.value);

        tdelayshock=xtrig/wactual*1000;

        tdelaymembrane=(xint)/wactual*1000+(xtrig-xint)/u2*1000;

        document.Timing.wactual.value = wactual;

        document.Timing.mactual.value = mactual;

        document.Timing.tdelayshock.value = tdelayshock;

        document.Timing.tdelaymembrane.value = tdelaymembrane;

    } catch (error) {
        console.error('Error in timing calculation:', error);
        alert('An error occurred during timing calculation: ' + error.message +
              '\n\nPlease ensure the shock calculation has been run first.');
    }
}



function calcBoost() {
    try {
        var drivervol,boostvol,prupture,percless,percmore,ru;

        var r,mw4,g4,psig,ppa,ppreboost,ppostboost,pboosttank;

        var deltap,deltam,t;

        psig=parseFloat(14.696);

        ppa=parseFloat(101325);

        drivervol = parseFloat(document.Boosttank.drivervol.value);

        boostvol = parseFloat(document.Boosttank.boostvol.value);

        prupture = parseFloat(document.Boosttank.prupture.value)*ppa/psig;

        percless = parseFloat(document.Boosttank.percless.value);

        percmore = parseFloat(document.Boosttank.percmore.value);

     

        // Get driver gas properties (boost calculations don't support custom mixtures)
        var driverGasId = document.Shock.drivergas.value;

        if (driverGasId === 'custom') {
            alert('Error: Boost tank calculations do not support custom gas mixtures. Please select a standard gas.');
            return;
        }

        var boostDriverGas = getGasById(driverGasId);
        if (!boostDriverGas) {
            alert('Error: Invalid driver gas selection for boost calculation');
            return;
        }

        g4 = boostDriverGas.gamma;
        mw4 = boostDriverGas.mw;

        ru=parseFloat(8314.510);

        // document.write("ru " + ru + "  ");

        r=ru/mw4;

        ppreboost=prupture*(1-percless/100);

        ppostboost=prupture*(1+percmore/100);

        deltap=ppostboost-ppreboost;

        t=parseFloat(document.Shock.tid.value);
        
        deltam=deltap*drivervol/r/t;


        pboosttank=deltam/boostvol*r*t+prupture;

       // document.write("pboosttank " + pboosttank + "  ");

       document.Boosttank.ppreboost.value = ppreboost*psig/ppa

       document.Boosttank.ppostboost.value = ppostboost*psig/ppa

       document.Boosttank.pboosttank.value = pboosttank*psig/ppa

    } catch (error) {
        console.error('Error in boost tank calculation:', error);
        alert('An error occurred during boost tank calculation: ' + error.message +
              '\n\nPlease check your input values and try again.');
    }
}

// Gas database with properties
var gasDatabase = {
    "1": {name: "Air", gamma: 1.4020, mw: 28.9660},
    "2": {name: "Argon", gamma: 1.6700, mw: 39.948000},
    "3": {name: "Carbon Dioxide", gamma: 1.29700, mw: 44.01000},
    "4": {name: "Helium", gamma: 1.66700, mw: 4.002602000},
    "5": {name: "Neon", gamma: 1.66700, mw: 20.179700},
    "6": {name: "Nitrogen", gamma: 1.40100, mw: 28.0134000},
    "7": {name: "Xenon", gamma: 1.66700, mw: 131.29000},
    "8": {name: "Sulfur Hexafluoride", gamma: 1.09200, mw: 146.06000},
    "10": {name: "Acetone Vapor", gamma: 1.14000, mw: 58.08000}
};

// Store mixture compositions
var mixtures = {
    driven: [],
    driver: []
};

// Custom gases stored in localStorage
var customGases = {};
var customGasCounter = 0;

// Flag to prevent saving state during initial load
var isLoadingState = false;

// Load custom gases from localStorage on page load
function loadCustomGases() {
    try {
        var stored = localStorage.getItem('wistl_custom_gases');
        if (stored) {
            var gases = JSON.parse(stored);
            customGases = {};
            for (var i = 0; i < gases.length; i++) {
                var gas = gases[i];
                customGases[gas.id] = gas;
                // Update counter to avoid ID collisions
                var numPart = parseInt(gas.id.replace('custom_', ''));
                if (numPart >= customGasCounter) {
                    customGasCounter = numPart + 1;
                }
            }
        }
    } catch (e) {
        console.error('Error loading custom gases:', e);
    }
}

// Save custom gases to localStorage
function saveCustomGases() {
    try {
        var gasArray = [];
        for (var id in customGases) {
            gasArray.push(customGases[id]);
        }
        localStorage.setItem('wistl_custom_gases', JSON.stringify(gasArray));
    } catch (e) {
        alert('Error saving custom gases. localStorage may be corrupted or disabled.\n\nYour custom gas will work in this session but may not persist after refresh.\n\nTo fix: Clear browser data (Settings > Privacy & Security > Cookies and Site Data > Clear Data)');
        console.error('Error saving custom gases:', e);
    }
}

// Get all gases (predefined + custom)
function getAllGases() {
    var allGases = {};
    // Add predefined gases
    for (var id in gasDatabase) {
        allGases[id] = gasDatabase[id];
    }
    // Add custom gases
    for (var id in customGases) {
        allGases[id] = customGases[id];
    }
    return allGases;
}

// Clear all output fields
function clearOutputFields() {
    // Shock calculator outputs
    document.Shock.pdriv.value = '';
    document.Shock.tdriv.value = '';
    document.Shock.ddriv.value = '';
    document.Shock.adriv.value = '';
    document.Shock.pshocked.value = '';
    document.Shock.tshocked.value = '';
    document.Shock.dshocked.value = '';
    document.Shock.ashocked.value = '';
    document.Shock.ushocked.value = '';
    document.Shock.pfour.value = '';
    document.Shock.pfourpsig.value = '';
    document.Shock.tfour.value = '';
    document.Shock.dfour.value = '';
    document.Shock.afour.value = '';
    document.Shock.preflect.value = '';
    document.Shock.treflect.value = '';
    document.Shock.dreflect.value = '';
    document.Shock.areflect.value = '';
    document.Shock.wreflect.value = '';
    document.Shock.mreflect.value = '';
    document.Shock.wshocked.value = '';
    document.Shock.outmw.value = '';
    document.Shock.outgamma.value = '';
    
    // Timing calculator outputs
    document.Timing.wactual.value = '';
    document.Timing.mactual.value = '';
    document.Timing.tdelayshock.value = '';
    document.Timing.tdelaymembrane.value = '';
    
    // Boost tank calculator outputs
    document.Boosttank.ppreboost.value = '';
    document.Boosttank.ppostboost.value = '';
    document.Boosttank.pboosttank.value = '';
}

// Check if localStorage is working
function checkLocalStorage() {
    try {
        var testKey = '__localStorage_test__';
        localStorage.setItem(testKey, 'test');
        localStorage.removeItem(testKey);
        return true;
    } catch (e) {
        return false;
    }
}

// Show localStorage error banner
function showLocalStorageWarning() {
    var banner = document.createElement('div');
    banner.style.cssText = 'background-color:#fff3cd;border:2px solid #ffc107;padding:15px;margin:10px;border-radius:5px;text-align:center;';
    banner.innerHTML = '<strong>⚠️ Storage Warning:</strong> Browser storage is corrupted or disabled. ' +
        'The calculator will work, but settings and custom gases won\'t be saved after refresh. ' +
        '<a href="LOCALSTORAGE_CORRUPTION_FIX.md" target="_blank" style="color:#0066cc;text-decoration:underline;">How to fix</a>';
    var container = document.getElementById('container') || document.body;
    container.insertBefore(banner, container.firstChild);
}

// Initialize custom gases on page load
window.onload = function() {
    loadCustomGases();
    updateCustomGasList();
    loadCalculatorState();
    clearOutputFields();

    // Check if localStorage is working and show warning if not
    if (!checkLocalStorage()) {
        showLocalStorageWarning();
    }

    // Add reset event listener to clear output fields
    document.Shock.addEventListener('reset', function() {
        // Use setTimeout to ensure the form reset completes first
        setTimeout(clearOutputFields, 0);
    });
};

// Save calculator state to localStorage
function saveCalculatorState() {
    try {
        var state = {
            mach: document.Shock.mach.value,
            pinit: document.Shock.pinit.value,
            tinit: document.Shock.tinit.value,
            tid: document.Shock.tid.value,
            areafact: document.Shock.areafact.value,
            drivengas: document.Shock.drivengas.value,
            drivergas: document.Shock.drivergas.value,
            // Save calculated outputs for import into x-t diagram
            pfour: document.Shock.pfour.value,
            pshocked: document.Shock.pshocked.value,
            ushocked: document.Shock.ushocked.value
        };
        
        // If custom mixture is selected, save the mixture compositions
        if (state.drivengas === 'custom') {
            state.drivenMixture = [];
            for (var i = 0; i < mixtures.driven.length; i++) {
                var index = mixtures.driven[i];
                var gasSelect = document.getElementById('driven_gas_' + index);
                var fracInput = document.getElementById('driven_frac_' + index);
                if (gasSelect && fracInput && gasSelect.value) {
                    state.drivenMixture.push({
                        gasId: gasSelect.value,
                        fraction: fracInput.value
                    });
                }
            }
        }
        
        if (state.drivergas === 'custom') {
            state.driverMixture = [];
            for (var i = 0; i < mixtures.driver.length; i++) {
                var index = mixtures.driver[i];
                var gasSelect = document.getElementById('driver_gas_' + index);
                var fracInput = document.getElementById('driver_frac_' + index);
                if (gasSelect && fracInput && gasSelect.value) {
                    state.driverMixture.push({
                        gasId: gasSelect.value,
                        fraction: fracInput.value
                    });
                }
            }
        }
        
        localStorage.setItem('wistl_calculator_state', JSON.stringify(state));
    } catch (e) {
        // Silently handle localStorage errors (might be corrupted, disabled, or full)
        // Don't spam the console on every keystroke if localStorage is broken
        if (!window._localStorageWarningShown) {
            console.warn('Cannot save calculator state - localStorage may be disabled or corrupted.');
            window._localStorageWarningShown = true;
        }
    }
}

// Load calculator state from localStorage
function loadCalculatorState() {
    try {
        var stored = localStorage.getItem('wistl_calculator_state');
        if (!stored) return;

        // Set flag to prevent saving during restoration
        isLoadingState = true;

        var state = JSON.parse(stored);
        
        // Restore basic inputs
        if (state.mach) document.Shock.mach.value = state.mach;
        if (state.pinit) document.Shock.pinit.value = state.pinit;
        if (state.tinit) document.Shock.tinit.value = state.tinit;
        if (state.tid) document.Shock.tid.value = state.tid;
        if (state.areafact) document.Shock.areafact.value = state.areafact;
        
        // Restore gas selections
        if (state.drivengas) {
            document.Shock.drivengas.value = state.drivengas;
            toggleMixtureInputs('driven');
            
            // Restore custom mixture if it was selected
            if (state.drivengas === 'custom' && state.drivenMixture && state.drivenMixture.length > 0) {
                // Clear any existing mixture components
                mixtures.driven = [];
                document.getElementById('driven_mixture_inputs').innerHTML = '';
                
                // Recreate mixture components
                for (var i = 0; i < state.drivenMixture.length; i++) {
                    addGasToMixture('driven');
                    var index = mixtures.driven[mixtures.driven.length - 1];
                    var gasSelect = document.getElementById('driven_gas_' + index);
                    var fracInput = document.getElementById('driven_frac_' + index);
                    
                    if (gasSelect && fracInput) {
                        gasSelect.value = state.drivenMixture[i].gasId;
                        fracInput.value = state.drivenMixture[i].fraction;
                    }
                }
                updateMixtureTotal('driven');
            }
        }
        
        if (state.drivergas) {
            document.Shock.drivergas.value = state.drivergas;
            toggleMixtureInputs('driver');
            
            // Restore custom mixture if it was selected
            if (state.drivergas === 'custom' && state.driverMixture && state.driverMixture.length > 0) {
                // Clear any existing mixture components
                mixtures.driver = [];
                document.getElementById('driver_mixture_inputs').innerHTML = '';
                
                // Recreate mixture components
                for (var i = 0; i < state.driverMixture.length; i++) {
                    addGasToMixture('driver');
                    var index = mixtures.driver[mixtures.driver.length - 1];
                    var gasSelect = document.getElementById('driver_gas_' + index);
                    var fracInput = document.getElementById('driver_frac_' + index);
                    
                    if (gasSelect && fracInput) {
                        gasSelect.value = state.driverMixture[i].gasId;
                        fracInput.value = state.driverMixture[i].fraction;
                    }
                }
                updateMixtureTotal('driver');
            }
        }
    } catch (e) {
        console.warn('Could not load calculator state (localStorage may be corrupted). Starting with default values.');
        // Try to clear corrupted data, but don't fail if this also errors
        try {
            localStorage.removeItem('wistl_calculator_state');
        } catch (clearError) {
            // If localStorage is completely corrupted at browser level, suggest manual fix
            console.warn('To fix localStorage corruption: Open Settings > Privacy & Security > Cookies and Site Data > Clear Data, or try Private Browsing mode.');
        }
    } finally {
        // Always clear the flag, even if there was an error
        isLoadingState = false;
    }
}

// Toggle Custom Gas Library visibility
function toggleCustomGasLibrary() {
    var section = document.getElementById('custom_gas_management');
    var button = document.getElementById('toggleCustomGasBtn');
    
    if (section.style.display === 'none' || section.style.display === '') {
        section.style.display = 'block';
        button.innerHTML = 'Manage Custom Gases ▲';
    } else {
        section.style.display = 'none';
        button.innerHTML = 'Manage Custom Gases ▼';
    }
}

// Show add gas form
function showAddGasForm() {
    var form = document.getElementById('add_gas_form');
    form.style.display = 'block';
    document.getElementById('new_gas_name').focus();
}

// Cancel add gas form
function cancelAddGasForm() {
    var form = document.getElementById('add_gas_form');
    form.style.display = 'none';
    document.getElementById('new_gas_name').value = '';
    document.getElementById('new_gas_gamma').value = '';
    document.getElementById('new_gas_mw').value = '';
}

// Add a new custom gas from the library
function addNewCustomGas() {
    var nameInput = document.getElementById('new_gas_name');
    var gammaInput = document.getElementById('new_gas_gamma');
    var mwInput = document.getElementById('new_gas_mw');
    
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
    
    // Update the list
    updateCustomGasList();
    
    // Refresh all dropdowns to include the new gas
    refreshAllGasDropdowns();
    
    // Clear and hide the form
    cancelAddGasForm();
    
    alert('Custom gas "' + name + '" saved successfully!');
}

// Show edit form for a gas
function editCustomGas(gasId) {
    var gas = customGases[gasId];
    if (!gas) return;
    
    var row = document.getElementById('gas_row_' + gasId);
    if (!row) return;
    
    // Store original values in case of cancel
    row.setAttribute('data-original-name', gas.name);
    row.setAttribute('data-original-gamma', gas.gamma);
    row.setAttribute('data-original-mw', gas.mw);
    
    // Replace row content with editable fields
    row.innerHTML = 
        '<td><input type="text" id="edit_name_' + gasId + '" value="' + gas.name + '" size="15"></td>' +
        '<td><input type="text" id="edit_gamma_' + gasId + '" value="' + gas.gamma + '" size="8"></td>' +
        '<td><input type="text" id="edit_mw_' + gasId + '" value="' + gas.mw + '" size="8"></td>' +
        '<td>' +
        '<button type="button" onclick="saveEditedCustomGas(\'' + gasId + '\')">Save</button> ' +
        '<button type="button" onclick="cancelEditCustomGas(\'' + gasId + '\')">Cancel</button>' +
        '</td>';
}

// Save edited gas properties
function saveEditedCustomGas(gasId) {
    var gas = customGases[gasId];
    if (!gas) return;
    
    var nameInput = document.getElementById('edit_name_' + gasId);
    var gammaInput = document.getElementById('edit_gamma_' + gasId);
    var mwInput = document.getElementById('edit_mw_' + gasId);
    
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
    
    // Check for duplicate name (excluding current gas)
    var allGases = getAllGases();
    for (var id in allGases) {
        if (id !== gasId && allGases[id].name.toLowerCase() === name.toLowerCase()) {
            alert('A gas with this name already exists. Please use a different name.');
            nameInput.focus();
            return;
        }
    }
    
    // Check if gas is currently in use and properties have changed significantly
    var propertiesChanged = Math.abs(gas.gamma - gamma) > 0.001 || Math.abs(gas.mw - mw) > 0.001;
    if (propertiesChanged) {
        var inUse = false;
        var types = ['driven', 'driver'];
        for (var t = 0; t < types.length; t++) {
            var type = types[t];
            for (var i = 0; i < mixtures[type].length; i++) {
                var index = mixtures[type][i];
                var select = document.getElementById(type + '_gas_' + index);
                if (select && select.value === gasId) {
                    inUse = true;
                    break;
                }
            }
            if (inUse) break;
        }
        
        if (inUse) {
            if (!confirm('Warning: This gas is currently selected in a mixture.\n\nChanging its properties will affect all calculations using it.\n\nDo you want to continue?')) {
                return;
            }
        }
    }
    
    // Update the gas properties
    gas.name = name;
    gas.gamma = gamma;
    gas.mw = mw;
    
    // Save to localStorage
    saveCustomGases();
    
    // Update the list display
    updateCustomGasList();
    
    // Update all dropdowns to reflect the new name
    refreshAllGasDropdowns();
    
    alert('Custom gas "' + name + '" updated successfully!');
}

// Cancel editing a gas
function cancelEditCustomGas(gasId) {
    // Simply refresh the list to restore original display
    updateCustomGasList();
}

// Update the custom gases list display
function updateCustomGasList() {
    var listDiv = document.getElementById('custom_gas_list');
    
    if (!listDiv) return;
    
    var gasCount = 0;
    for (var id in customGases) {
        gasCount++;
    }
    
    if (gasCount === 0) {
        listDiv.innerHTML = '<span style="color:#6c757d; font-style:italic;">No custom gases defined yet.</span>';
        return;
    }
    
    var html = '<table>';
    html += '<tr>';
    html += '<td>Gas Name</td>';
    html += '<td>Gamma (γ)</td>';
    html += '<td>MW (g/mol)</td>';
    html += '<td>Actions</td>';
    html += '</tr>';
    
    for (var id in customGases) {
        var gas = customGases[id];
        html += '<tr id="gas_row_' + id + '">';
        html += '<td>★ ' + gas.name + '</td>';
        html += '<td>' + gas.gamma.toFixed(4) + '</td>';
        html += '<td>' + gas.mw.toFixed(4) + '</td>';
        html += '<td>';
        html += '<button type="button" onclick="editCustomGas(\'' + id + '\')">Edit</button> ';
        html += '<button type="button" onclick="deleteCustomGas(\'' + id + '\')">Delete</button>';
        html += '</td>';
        html += '</tr>';
    }
    
    html += '</table>';
    listDiv.innerHTML = html;
}

// Delete a custom gas
function deleteCustomGas(gasId) {
    var gas = customGases[gasId];
    if (!gas) return;
    
    if (!confirm('Are you sure you want to delete "' + gas.name + '"?\n\nThis will remove it from all mixture dropdowns.')) {
        return;
    }
    
    // Check if this gas is currently being used in any mixture
    var inUse = false;
    var types = ['driven', 'driver'];
    for (var t = 0; t < types.length; t++) {
        var type = types[t];
        for (var i = 0; i < mixtures[type].length; i++) {
            var index = mixtures[type][i];
            var select = document.getElementById(type + '_gas_' + index);
            if (select && select.value === gasId) {
                inUse = true;
                break;
            }
        }
        if (inUse) break;
    }
    
    if (inUse) {
        alert('Warning: This gas is currently selected in a mixture. It will be removed from the selection.');
    }
    
    // Delete the gas
    delete customGases[gasId];
    saveCustomGases();
    
    // Update the list
    updateCustomGasList();
    
    // Remove from all dropdowns and reset selections if this gas was selected
    for (var t = 0; t < types.length; t++) {
        var type = types[t];
        for (var i = 0; i < mixtures[type].length; i++) {
            var index = mixtures[type][i];
            var select = document.getElementById(type + '_gas_' + index);
            if (select) {
                // If this gas was selected, reset the selection
                if (select.value === gasId) {
                    select.value = '';
                }
                
                // Remove the option from the dropdown
                var options = select.options;
                for (var j = options.length - 1; j >= 0; j--) {
                    if (options[j].value === gasId) {
                        select.remove(j);
                        break;
                    }
                }
            }
        }
    }
    
    // Update mixture totals in case we removed a selected gas
    updateMixtureTotal('driven');
    updateMixtureTotal('driver');
}

// Toggle visibility of mixture input sections
function toggleMixtureInputs(type) {
    var selectElement = document.Shock[type + 'gas'];
    var divElement = document.getElementById(type + '_mixture_div');

    if (selectElement.value === 'custom') {
        divElement.style.display = 'block';
        // Initialize with one gas input if empty
        if (mixtures[type].length === 0) {
            addGasToMixture(type);
        }
    } else {
        divElement.style.display = 'none';
    }

    // Save state whenever gas selection changes (unless we're loading state)
    if (!isLoadingState) {
        saveCalculatorState();
    }
}

// Add a new gas component to the mixture
function addGasToMixture(type) {
    var container = document.getElementById(type + '_mixture_inputs');
    var index = mixtures[type].length;
    
    var div = document.createElement('div');
    div.id = type + '_mix_' + index;
    div.className = 'mixture-component';
    
    var select = document.createElement('select');
    select.id = type + '_gas_' + index;
    select.onchange = function() { 
        handleGasSelection(type, index);
    };
    
    // Add gas options
    var emptyOption = document.createElement('option');
    emptyOption.value = '';
    emptyOption.text = 'Select Gas';
    select.appendChild(emptyOption);
    
    // Add "Define New Gas..." option
    var defineOption = document.createElement('option');
    defineOption.value = 'define_new';
    defineOption.text = '--- Define New Gas... ---';
    select.appendChild(defineOption);
    
    // Collect all gases and sort alphabetically
    var allGasOptions = [];
    
    // Add predefined gases
    for (var gasId in gasDatabase) {
        allGasOptions.push({
            value: gasId,
            text: gasDatabase[gasId].name,
            sortKey: gasDatabase[gasId].name
        });
    }
    
    // Add custom gases (with special marker)
    for (var gasId in customGases) {
        allGasOptions.push({
            value: gasId,
            text: '★ ' + customGases[gasId].name,
            sortKey: customGases[gasId].name
        });
    }
    
    // Sort alphabetically by name
    allGasOptions.sort(function(a, b) {
        return a.sortKey.localeCompare(b.sortKey);
    });
    
    // Add sorted options to select
    for (var i = 0; i < allGasOptions.length; i++) {
        var option = document.createElement('option');
        option.value = allGasOptions[i].value;
        option.text = allGasOptions[i].text;
        select.appendChild(option);
    }
    
    var label = document.createElement('span');
    label.textContent = ' Mole Fraction: ';
    
    var input = document.createElement('input');
    input.type = 'text';
    input.id = type + '_frac_' + index;
    input.size = 6;
    input.maxLength = 8;
    input.placeholder = '0.0';
    input.onkeyup = function() { updateMixtureTotal(type); };
    input.onchange = function() { updateMixtureTotal(type); };
    
    var removeBtn = document.createElement('button');
    removeBtn.type = 'button';
    removeBtn.innerHTML = 'Remove';
    removeBtn.onclick = function() { removeGasFromMixture(type, index); };
    
    div.appendChild(select);
    div.appendChild(label);
    div.appendChild(input);
    div.appendChild(document.createTextNode(' '));
    div.appendChild(removeBtn);
    
    // Container for custom gas definition form (initially hidden)
    var defineDiv = document.createElement('div');
    defineDiv.id = type + '_define_' + index;
    defineDiv.style.display = 'none';
    defineDiv.style.marginLeft = '20px';
    defineDiv.style.marginTop = '8px';
    defineDiv.style.padding = '10px';
    defineDiv.style.backgroundColor = '#fff3cd';
    defineDiv.style.border = '1px solid #ffc107';
    defineDiv.style.borderRadius = '3px';
    div.appendChild(document.createElement('br'));
    div.appendChild(defineDiv);
    
    container.appendChild(div);
    mixtures[type].push(index);
    
    updateMixtureTotal(type);
}

// Handle gas selection (including "Define New Gas" option)
function handleGasSelection(type, index) {
    var select = document.getElementById(type + '_gas_' + index);
    var defineDiv = document.getElementById(type + '_define_' + index);

    if (select.value === 'define_new') {
        // Show the custom gas definition form
        showDefineGasForm(type, index);
    } else {
        // Hide the definition form if it was showing
        if (defineDiv) {
            defineDiv.style.display = 'none';
        }
        updateMixtureTotal(type);
        // Save state when gas selection changes (unless we're loading state)
        if (!isLoadingState) {
            saveCalculatorState();
        }
    }
}

// Show the custom gas definition form
function showDefineGasForm(type, index) {
    var defineDiv = document.getElementById(type + '_define_' + index);
    
    defineDiv.innerHTML = '<b>Define New Gas:</b><br>' +
        'Name: <input type="text" id="' + type + '_newgas_name_' + index + '" placeholder="e.g., Methane" style="width:140px;"> ' +
        'Gamma (γ): <input type="text" id="' + type + '_newgas_gamma_' + index + '" placeholder="e.g., 1.31" style="width:80px;"> ' +
        'MW (g/mol): <input type="text" id="' + type + '_newgas_mw_' + index + '" placeholder="e.g., 16.04" style="width:80px;"> ' +
        '<button type="button" onclick="saveCustomGas(\'' + type + '\', ' + index + ')">Save Gas</button> ' +
        '<button type="button" onclick="cancelDefineGas(\'' + type + '\', ' + index + ')">Cancel</button>';
    
    defineDiv.style.display = 'block';
}

// Save a new custom gas
function saveCustomGas(type, index) {
    var nameInput = document.getElementById(type + '_newgas_name_' + index);
    var gammaInput = document.getElementById(type + '_newgas_gamma_' + index);
    var mwInput = document.getElementById(type + '_newgas_mw_' + index);
    
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
    
    // Update the gas selection dropdown
    var select = document.getElementById(type + '_gas_' + index);
    var option = document.createElement('option');
    option.value = gasId;
    option.text = '★ ' + name;
    option.selected = true;
    select.appendChild(option);
    
    // Hide the definition form
    var defineDiv = document.getElementById(type + '_define_' + index);
    defineDiv.style.display = 'none';
    
    // Update custom gases list
    updateCustomGasList();
    
    // Refresh all other dropdowns to include the new gas
    refreshAllGasDropdowns();
    
    alert('Custom gas "' + name + '" saved successfully!');
    updateMixtureTotal(type);
}

// Cancel custom gas definition
function cancelDefineGas(type, index) {
    var select = document.getElementById(type + '_gas_' + index);
    var defineDiv = document.getElementById(type + '_define_' + index);
    
    select.value = '';
    defineDiv.style.display = 'none';
    updateMixtureTotal(type);
}

// Refresh all gas dropdowns to include newly added custom gases
function refreshAllGasDropdowns() {
    var types = ['driven', 'driver'];
    
    for (var t = 0; t < types.length; t++) {
        var type = types[t];
        for (var i = 0; i < mixtures[type].length; i++) {
            var index = mixtures[type][i];
            var select = document.getElementById(type + '_gas_' + index);
            if (select) {
                var currentValue = select.value;
                
                // Remove old custom gas options
                var options = select.options;
                for (var j = options.length - 1; j >= 0; j--) {
                    if (options[j].value.indexOf('custom_') === 0) {
                        select.remove(j);
                    }
                }
                
                // Add updated custom gas options
                for (var gasId in customGases) {
                    var option = document.createElement('option');
                    option.value = gasId;
                    option.text = '★ ' + customGases[gasId].name;
                    select.appendChild(option);
                }
                
                // Restore previous selection if it still exists
                select.value = currentValue;
            }
        }
    }
}

// Remove a gas component from the mixture
function removeGasFromMixture(type, index) {
    var div = document.getElementById(type + '_mix_' + index);
    if (div) {
        div.remove();
        var idx = mixtures[type].indexOf(index);
        if (idx > -1) {
            mixtures[type].splice(idx, 1);
        }
        updateMixtureTotal(type);
        // Save state when component is removed (unless we're loading state)
        if (!isLoadingState) {
            saveCalculatorState();
        }
    }
}

// Update and display the total mole fraction
function updateMixtureTotal(type) {
    var total = 0.0;
    var valid = true;

    for (var i = 0; i < mixtures[type].length; i++) {
        var index = mixtures[type][i];
        var input = document.getElementById(type + '_frac_' + index);
        if (input) {
            var val = parseFloat(input.value);
            if (!isNaN(val)) {
                total += val;
            }
        }
    }

    var totalSpan = document.getElementById(type + '_total');
    var displayTotal = Math.round(total * 100 * 100) / 100; // Round to 2 decimal places

    if (Math.abs(total - 1.0) < 0.001) {
        totalSpan.innerHTML = '<span style="color:#28a745;">Total: ' + displayTotal + '% ✓</span>';
    } else if (total > 0) {
        totalSpan.innerHTML = '<span style="color:#dc3545;">Total: ' + displayTotal + '% (must equal 100%)</span>';
    } else {
        totalSpan.innerHTML = '<span style="color:#6c757d;">Total: 0%</span>';
    }

    // Save state when fractions change (unless we're loading state)
    if (!isLoadingState) {
        saveCalculatorState();
    }
}

// Calculate mixture properties (gamma and molecular weight)
function calculateMixtureProperties(type) {
    var totalFraction = 0.0;
    var components = [];
    
    // Collect all components
    for (var i = 0; i < mixtures[type].length; i++) {
        var index = mixtures[type][i];
        var gasSelect = document.getElementById(type + '_gas_' + index);
        var fracInput = document.getElementById(type + '_frac_' + index);
        
        if (gasSelect && fracInput) {
            var gasId = gasSelect.value;
            var fraction = parseFloat(fracInput.value);
            
            if (gasId && !isNaN(fraction) && fraction > 0) {
                components.push({
                    gasId: gasId,
                    fraction: fraction
                });
                totalFraction += fraction;
            }
        }
    }
    
    // Validate total fraction
    if (Math.abs(totalFraction - 1.0) > 0.001) {
        alert('Error: Mole fractions must sum to 1.0 (100%). Current total: ' + (totalFraction * 100).toFixed(2) + '%');
        return null;
    }
    
    if (components.length === 0) {
        alert('Error: No gas components specified in custom mixture.');
        return null;
    }
    
    // Get all available gases (predefined + custom)
    var allGases = getAllGases();
    
    // Calculate mixture molecular weight (simple weighted average)
    var mw_mix = 0.0;
    for (var i = 0; i < components.length; i++) {
        var gas = allGases[components[i].gasId];
        if (!gas) continue;
        mw_mix += components[i].fraction * gas.mw;
    }
    
    // Calculate mixture gamma using specific heats
    // For each gas: Cv = R/(gamma - 1), Cp = gamma * Cv
    // Then: Cv_mix = sum(x_i * Cv_i), Cp_mix = sum(x_i * Cp_i)
    // Finally: gamma_mix = Cp_mix / Cv_mix
    
    var cv_mix = 0.0;
    var cp_mix = 0.0;
    
    for (var i = 0; i < components.length; i++) {
        var gas = allGases[components[i].gasId];
        if (!gas) continue;
        var x = components[i].fraction;
        var gamma = gas.gamma;
        
        // Specific heats (normalized by R/MW, which cancels out in the ratio)
        var cv = 1.0 / (gamma - 1.0);
        var cp = gamma * cv;
        
        cv_mix += x * cv;
        cp_mix += x * cp;
    }
    
    var gamma_mix = cp_mix / cv_mix;
    
    return {
        gamma: gamma_mix,
        mw: mw_mix
    };
}


function showEquationsModal() {
    document.getElementById('equationsModal').style.display = 'block';
    // Trigger MathJax to render the equations
    if (window.MathJax) {
        MathJax.typesetPromise();
    }
}

function closeEquationsModal() {
    document.getElementById('equationsModal').style.display = 'none';
}

// Close modal when clicking outside of it
window.onclick = function(event) {
    var modal = document.getElementById('equationsModal');
    if (event.target == modal) {
        closeEquationsModal();
    }
}
