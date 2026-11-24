/**
 * Shock tube, timing and boost tank calculators, the custom gas library
 * panel, and persistence of the shock tube inputs (browser storage and
 * shareable links).
 */

const STATE_STORAGE_KEY = 'wistl_calculator_state';
const PSI_PER_PA = 14.696 / 101325;

/** Results other tools read. Never read values back out of output fields. */
const appState = {
    shock: null,        // result of computeShockTube for the current inputs
    shockInputs: null,  // {driven, driver, p1, t1, t4} that produced it
    shockListeners: []  // called whenever appState.shock changes
};

let drivenPicker, driverPicker;

/* ------------------------------------------------------------------ */
/* Shock tube                                                          */
/* ------------------------------------------------------------------ */

const SHOCK_OUTPUT_IDS = [
    'pdriv', 'pshocked', 'pfour', 'preflect', 'tdriv', 'tshocked', 'tfour', 'treflect',
    'ddriv', 'dshocked', 'dfour', 'dreflect', 'adriv', 'ashocked', 'afour', 'areflect',
    'gdriv', 'gshocked', 'gfour', 'greflect',
    'ushocked', 'outmach', 'wshocked', 'wreflect', 'mreflect', 'pfourpsig', 'outmw', 'outtailor', 'outtailormach'
];

function shockMode() {
    return byId('shockForm').elements.mode.value;
}

function updateShockMode() {
    const mode = shockMode();
    byId('shockForm').querySelectorAll('[data-mode]').forEach(function (node) {
        node.hidden = node.dataset.mode !== mode;
    });
}

function setShockResult(result, inputs) {
    appState.shock = result;
    appState.shockInputs = result ? inputs : null;
    appState.shockListeners.forEach(function (listener) { listener(); });
}

function clearShockOutputs() {
    SHOCK_OUTPUT_IDS.forEach(function (id) { byId(id).value = ''; });
    byId('shockOutput').classList.remove('stale');
    setShockResult(null);
}

/** Write a result to the output fields in the current display units. */
function displayShockResult(result) {
    const s1 = result.state1, s2 = result.state2, s4 = result.state4, s5 = result.state5;
    const set = function (id, text) { byId(id).value = text; };
    const column = function (suffix, s) {
        set('p' + suffix, formatQuantity('pressure', s.p));
        set('t' + suffix, formatQuantity('temperature', s.T));
        set('d' + suffix, formatSig(s.rho));
        set('a' + suffix, formatSig(s.a));
        set('g' + suffix, formatSig(s.gamma, 5));
    };
    column('driv', s1);
    column('shocked', s2);
    column('four', s4);
    column('reflect', s5);
    set('ushocked', formatSig(s2.u));
    set('outmach', formatSig(result.mach));
    set('wshocked', formatSig(result.shockSpeed));
    set('wreflect', formatSig(result.reflectedShockSpeed));
    set('mreflect', formatSig(result.reflectedMach));
    // The legacy estimate is stored in psi gauge
    set('pfourpsig', result.driverPressurePsig === null ? '—'
        : formatQuantity('gauge', result.driverPressurePsig / PSI_PER_PA));
    set('outmw', formatSig(result.mw));
    set('outtailor', formatSig(result.tailoringRatio, 4));
    set('outtailormach', formatSig(result.tailoredMach, 4));
    byId('outmwLabel').title = result.pinf > 0
        ? 'Effective value for the stiffened gas; p∞ = ' + formatSig(result.pinf / 1e6, 4) + ' MPa'
        : '';
}

function calcShock() {
    const form = byId('shockForm');
    const messages = byId('shockMessages');
    clearInvalid(form);

    let result, inputs;
    try {
        inputs = {
            driven: drivenPicker.resolve(),
            driver: driverPicker.resolve(),
            p1: readQuantity(byId('pinit'), 'Driven pressure'),
            t1: readQuantity(byId('tinit'), 'Driven temperature'),
            t4: readQuantity(byId('tid'), 'Driver temperature'),
            areaRatio: readNumber(byId('arearatio'), 'Area ratio', { min: 1 }),
            areaFactor: readNumber(byId('areafact'), 'Empirical area factor', { positive: true }),
            variableGamma: byId('gasModel').value === 'nasa',
            formatTemperature: function (T) { return formatQuantity('temperature', T, 4) + ' ' + Units.label('temperature'); }
        };

        let mach;
        if (shockMode() === 'pressure') {
            const p4 = readQuantity(byId('p4input'), 'Driver pressure');
            if (p4 <= inputs.p1) throw new FieldError(byId('p4input'), 'Driver pressure must exceed the driven pressure to produce a shock.');
            mach = machFromDriverPressure({ ...inputs, p4: p4 });
        } else {
            mach = readNumber(byId('mach'), 'Mach number');
            if (mach <= 1) throw new FieldError(byId('mach'), 'Mach number must be greater than 1 for a shock wave.');
        }
        result = computeShockTube({ ...inputs, mach: mach });
        result.tailoredMach = tailoredMach(inputs);
    } catch (error) {
        // Never leave results on screen that do not belong to the inputs
        clearShockOutputs();
        reportError(messages, error);
        return;
    }

    displayShockResult(result);
    showMessages(messages, result.warnings.map(function (w) { return { type: 'warning', text: w.message }; }));
    byId('shockOutput').classList.remove('stale');
    setShockResult(result, inputs);
    saveShockState();
}

/* ------------------------------------------------------------------ */
/* Persistence: browser storage and shareable links                    */
/* ------------------------------------------------------------------ */

// Saved states and links hold SI values (Pa, K) whatever units are displayed
const SHOCK_FIELDS = { m: 'mach', p4: 'p4input', p1: 'pinit', t1: 'tinit', t4: 'tid', ar: 'arearatio', af: 'areafact' };

function readShockInputState() {
    const state = {
        mode: shockMode(), gm: byId('gasModel').value,
        driven: drivenPicker.getValue(), driver: driverPicker.getValue()
    };
    for (const key in SHOCK_FIELDS) {
        const input = byId(SHOCK_FIELDS[key]);
        const si = input.dataset.quantity ? Units.inputSI(input) : null;
        state[key] = si === null ? input.value.trim() : String(Number(si.toPrecision(12)));
    }
    return state;
}

function applyShockInputState(state) {
    if (state.mode === 'pressure' || state.mode === 'mach') byId('shockForm').elements.mode.value = state.mode;
    if (state.gm === 'nasa' || state.gm === 'constant') byId('gasModel').value = state.gm;
    for (const key in SHOCK_FIELDS) {
        if (state[key] === undefined || state[key] === null) continue;
        const input = byId(SHOCK_FIELDS[key]);
        const numeric = String(state[key]).trim() !== '' && isFinite(Number(state[key]));
        if (input.dataset.quantity && numeric) Units.setInput(input, Number(state[key]));
        else input.value = state[key];
    }
    if (state.driven) drivenPicker.setValue(state.driven);
    if (state.driver) driverPicker.setValue(state.driver);
    updateShockMode();
}

function shockStateToQuery(state) {
    const params = new URLSearchParams();
    params.set('mode', state.mode);
    if (state.gm === 'nasa') params.set('gm', 'nasa');
    for (const key in SHOCK_FIELDS) {
        const usedInMode = (key !== 'm' || state.mode === 'mach') && (key !== 'p4' || state.mode === 'pressure');
        if (usedInMode) params.set(key, state[key]);
    }
    params.set('dn', encodeGasValue(state.driven));
    params.set('dr', encodeGasValue(state.driver));
    return params.toString();
}

function shockStateFromQuery(query) {
    const params = new URLSearchParams(query);
    if (!params.has('dn') || !params.has('dr')) return null;
    const state = { mode: params.get('mode') || 'mach', gm: params.get('gm') || 'constant' };
    for (const key in SHOCK_FIELDS) {
        if (params.has(key)) state[key] = params.get(key);
    }
    state.driven = decodeGasValue(params.get('dn'));
    state.driver = decodeGasValue(params.get('dr'));
    return state;
}

/** Convert a state saved by an earlier version of the page. */
function migrateStoredState(stored) {
    if (stored.driven || !stored.drivengas) return stored;
    return {
        mode: 'mach',
        m: stored.mach, p1: stored.pinit, t1: stored.tinit, t4: stored.tid, af: stored.areafact,
        driven: { gasId: stored.drivengas, mixture: stored.drivenMixture || [] },
        driver: { gasId: stored.drivergas, mixture: stored.driverMixture || [] }
    };
}

function saveShockState() {
    const state = readShockInputState();
    try {
        localStorage.setItem(STATE_STORAGE_KEY, JSON.stringify(state));
    } catch (e) { /* storage unavailable: the page still works */ }
    history.replaceState(null, '', '?' + shockStateToQuery(state) + location.hash);
}

/** Links refer to NASA database species by name ("n:C3H8"). */
function linkNeedsNasaDatabase() {
    const params = new URLSearchParams(location.search);
    return /(^|,)n:/.test((params.get('dn') || '') + ',' + (params.get('dr') || ''));
}

/** @returns {boolean} true when the state came from a shared link */
function loadShockState() {
    try {
        const fromLink = location.search.length > 1 ? shockStateFromQuery(location.search.slice(1)) : null;
        if (fromLink) {
            applyShockInputState(fromLink);
            return true;
        }
    } catch (e) {
        console.warn('Ignoring unreadable link parameters.', e);
    }
    try {
        const stored = JSON.parse(localStorage.getItem(STATE_STORAGE_KEY) || 'null');
        if (stored) applyShockInputState(migrateStoredState(stored));
    } catch (e) {
        console.warn('Ignoring unreadable saved state.', e);
    }
    return false;
}

function copyShockLink() {
    const url = location.origin + location.pathname + '?' + shockStateToQuery(readShockInputState());
    const done = function () {
        showMessages(byId('shockMessages'), [{ type: 'success', text: 'Link copied.' }]);
    };
    const fallback = function () {
        history.replaceState(null, '', url);
        showMessages(byId('shockMessages'), [{ type: 'info', text: 'Link is in the address bar.' }]);
    };
    if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(url).then(done, fallback);
    } else {
        fallback();
    }
}

/* ------------------------------------------------------------------ */
/* Timing                                                              */
/* ------------------------------------------------------------------ */

const TIMING_OUTPUT_IDS = ['wactual', 'mactual', 'tdelayshock', 'tdelaymembrane'];

function calcTime() {
    const form = byId('timingForm');
    const messages = byId('timingMessages');
    clearInvalid(form);
    TIMING_OUTPUT_IDS.forEach(function (id) { byId(id).value = ''; });

    try {
        const xwave = readNumber(byId('xwave'), 'Transducer spacing', { positive: true });
        const twave = readNumber(byId('twave'), 'Transit time', { positive: true }) / 1000;
        const xtrig = readNumber(byId('xtrig'), 'Trigger-to-window distance', { min: 0 });
        const xint = readNumber(byId('xint'), 'Trigger-to-interface distance', { min: 0 });

        const wactual = xwave / twave;
        byId('wactual').value = formatSig(wactual);
        byId('tdelayshock').value = formatSig(xtrig / wactual * 1000);

        const shock = appState.shock;
        if (!shock) {
            byId('mactual').value = '—';
            byId('tdelaymembrane').value = '—';
            showMessages(messages, [{
                type: 'warning',
                text: 'Run the shock tube calculation to get the Mach number and interface delay.'
            }]);
            return;
        }

        const a1 = shock.state1.a, u2 = shock.state2.u;
        byId('mactual').value = formatSig(wactual / a1);
        // Shock travels to the interface, then the interface moves with the shocked gas
        byId('tdelaymembrane').value = formatSig((xint / wactual + (xtrig - xint) / u2) * 1000);
        showMessages(messages, []);
    } catch (error) {
        reportError(messages, error);
    }
}

/* ------------------------------------------------------------------ */
/* Boost tank                                                          */
/* ------------------------------------------------------------------ */

const BOOST_OUTPUT_IDS = ['ppreboost', 'ppostboost', 'pboosttank'];

function calcBoost() {
    const form = byId('boostForm');
    const messages = byId('boostMessages');
    clearInvalid(form);
    BOOST_OUTPUT_IDS.forEach(function (id) { byId(id).value = ''; });

    try {
        const driverVolume = readNumber(byId('drivervol'), 'Driver volume', { positive: true });
        const boostVolume = readNumber(byId('boostvol'), 'Boost tank volume', { positive: true });
        // Gauge pressures in the display unit; the relations are linear in them
        const rupture = readNumber(byId('prupture'), 'Rupture pressure', { positive: true });
        const less = readNumber(byId('percless'), 'Pre-boost margin', { min: 0, max: 100 });
        const more = readNumber(byId('percmore'), 'Post-boost margin', { min: 0 });

        // Isothermal transfer: the boost tank ends at the rupture pressure after
        // raising the driver from the pre-boost to the post-boost pressure.
        const preBoost = rupture * (1 - less / 100);
        const postBoost = rupture * (1 + more / 100);
        const boostTank = (postBoost - preBoost) * driverVolume / boostVolume + rupture;

        byId('ppreboost').value = formatSig(preBoost);
        byId('ppostboost').value = formatSig(postBoost);
        byId('pboosttank').value = formatSig(boostTank);
        showMessages(messages, []);
    } catch (error) {
        reportError(messages, error);
    }
}

/* ------------------------------------------------------------------ */
/* Custom gas library panel                                            */
/* ------------------------------------------------------------------ */

function updateGasFormEos() {
    const eos = byId('gasFormEos').value;
    byId('gasForm').querySelectorAll('[data-eos]').forEach(function (node) {
        node.hidden = node.dataset.eos !== eos;
    });
}

function resetGasForm() {
    byId('gasForm').reset();
    byId('gasFormId').value = '';
    byId('gasFormSave').textContent = 'Add';
    byId('gasFormCancel').hidden = true;
    clearInvalid(byId('gasForm'));
    updateGasFormEos();
}

function editGasInForm(id) {
    const gas = GasLibrary.custom[id];
    if (!gas) return;
    byId('gasFormId').value = id;
    byId('gasFormName').value = gas.name;
    byId('gasFormGamma').value = gas.gamma;
    byId('gasFormEos').value = gas.pinf > 0 ? 'stiffened' : 'ideal';
    byId('gasFormMw').value = gas.pinf > 0 ? '' : gas.mw;
    byId('gasFormPinf').value = gas.pinf > 0 ? gas.pinf / 1e6 : '';
    byId('gasFormRho').value = gas.pinf > 0 ? formatSig(eosDensity(gas, STP_REF_P, STP_REF_T), 8) : '';
    byId('gasFormSave').textContent = 'Save changes';
    byId('gasFormCancel').hidden = false;
    updateGasFormEos();
    byId('gasFormName').focus();
}

function submitGasForm() {
    const messages = byId('gasMessages');
    clearInvalid(byId('gasForm'));
    try {
        const name = byId('gasFormName').value.trim();
        if (!name) throw new FieldError(byId('gasFormName'), 'Enter a name.');
        const gamma = readNumber(byId('gasFormGamma'), 'γ');
        if (gamma <= 1) throw new FieldError(byId('gasFormGamma'), 'γ must be greater than 1.');

        let mw, pinf = 0;
        if (byId('gasFormEos').value === 'stiffened') {
            pinf = readNumber(byId('gasFormPinf'), 'p∞', { positive: true }) * 1e6;
            const rho = readNumber(byId('gasFormRho'), 'Reference density', { positive: true });
            mw = stiffenedEffectiveMw(pinf, rho);
        } else {
            mw = readNumber(byId('gasFormMw'), 'Molecular weight', { positive: true });
        }

        const id = byId('gasFormId').value || undefined;
        GasLibrary.put({ name: name, gamma: gamma, mw: mw, pinf: pinf }, id);
        resetGasForm();
        const items = [{ type: 'success', text: (id ? 'Updated ' : 'Added ') + name + '.' }];
        if (!GasLibrary.storageWorks) {
            items.push({ type: 'warning', text: 'Browser storage is unavailable; custom gases will not be kept.' });
        }
        showMessages(messages, items);
    } catch (error) {
        reportError(messages, error);
    }
}

function gasInUse(id) {
    let used = false;
    GasPicker.instances.forEach(function (picker) {
        const value = picker.getValue();
        if (value.gasId === id) used = true;
        if (value.gasId === MIXTURE_VALUE && value.mixture.some(function (row) { return row.gasId === id; })) used = true;
    });
    return used;
}

function renderGasList() {
    const gases = Object.values(GasLibrary.custom);
    byId('gasList').replaceChildren(...gases.map(function (gas) {
        const props = gas.fit
            ? 'NASA fit, γ ' + formatSig(gas.gamma, 5) + ' at 298 K, ' + formatSig(gas.mw, 6) + ' g/mol'
            : gas.pinf > 0
            ? 'γ ' + formatSig(gas.gamma, 5) + ', p∞ ' + formatSig(gas.pinf / 1e6, 5) + ' MPa, ' +
              formatSig(eosDensity(gas, STP_REF_P, STP_REF_T), 5) + ' kg/m³'
            : 'γ ' + formatSig(gas.gamma, 5) + ', ' + formatSig(gas.mw, 6) + ' g/mol';
        return el('div', { class: 'gas-entry' },
            el('strong', {}, gas.name),
            el('span', { class: 'props' }, props),
            gas.fit ? null : el('button', { type: 'button', class: 'secondary small', onclick: function () { editGasInForm(gas.id); } }, 'Edit'),
            el('button', {
                type: 'button', class: 'danger small',
                onclick: function () {
                    const note = gasInUse(gas.id) ? ' It is currently selected and will be deselected.' : '';
                    if (!confirm('Delete "' + gas.name + '"?' + note)) return;
                    GasLibrary.remove(gas.id);
                    if (byId('gasFormId').value === gas.id) resetGasForm();
                }
            }, 'Delete'));
    }));
}

/* ------------------------------------------------------------------ */
/* NASA species search                                                 */
/* ------------------------------------------------------------------ */

const NASA_RESULT_LIMIT = 24;

function renderNasaResults(database) {
    const container = byId('nasaResults');
    const query = byId('nasaSearch').value;
    const matches = searchNasaSpecies(database, query);
    if (!query.trim()) {
        container.replaceChildren();
        return;
    }
    if (matches.length === 0) {
        container.replaceChildren('No species match. Try a formula such as C3H8.');
        return;
    }
    const chips = matches.slice(0, NASA_RESULT_LIMIT).map(function (species) {
        const added = GasLibrary.findSpecies(species) !== null;
        return el('button', {
            type: 'button', class: 'secondary small', disabled: added,
            title: formatSig(database[species].mw, 5) + ' g/mol' + (added ? ' (already added)' : ''),
            onclick: function () {
                GasLibrary.addSpecies(species, database[species]);
                showMessages(byId('gasMessages'), [{ type: 'success', text: 'Added ' + species + '.' }]);
                renderNasaResults(database);
            }
        }, species);
    });
    const more = matches.length - chips.length;
    container.replaceChildren(...chips, more > 0 ? '+' + more + ' more' : '');
}

function searchNasa() {
    loadNasaDatabase().then(renderNasaResults, function (error) {
        byId('nasaResults').replaceChildren(error.message);
    });
}

/* ------------------------------------------------------------------ */
/* Equations dialog                                                    */
/* ------------------------------------------------------------------ */

const MATHJAX_URL = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js';
let equationsOpener = null;

function showEquationsModal() {
    const modal = byId('equationsModal');
    equationsOpener = document.activeElement;
    modal.classList.add('open');
    byId('closeEquationsBtn').focus();
    if (window.MathJax && MathJax.typesetPromise) {
        MathJax.typesetPromise([modal]);
    } else if (!byId('MathJax-script')) {
        // Loaded on first use; MathJax typesets the page when it starts
        document.head.append(el('script', { id: 'MathJax-script', src: MATHJAX_URL, async: true }));
    }
}

function closeEquationsModal() {
    byId('equationsModal').classList.remove('open');
    if (equationsOpener && equationsOpener.focus) equationsOpener.focus();
}

/* ------------------------------------------------------------------ */
/* Start-up                                                            */
/* ------------------------------------------------------------------ */

function onSubmit(formId, handler) {
    byId(formId).addEventListener('submit', function (event) {
        event.preventDefault();
        handler();
    });
}

function onReset(formId, outputIds, messagesId) {
    byId(formId).addEventListener('reset', function () {
        clearInvalid(byId(formId));
        outputIds.forEach(function (id) { byId(id).value = ''; });
        showMessages(byId(messagesId), []);
    });
}

function setupUnitPickers() {
    const pickers = { pressure: byId('unitPressure'), temperature: byId('unitTemperature') };
    for (const quantity in pickers) {
        const select = pickers[quantity];
        select.replaceChildren(...Object.keys(Units.table[quantity]).map(function (unit) {
            return el('option', { value: unit }, unit);
        }));
        select.value = Units.current[quantity];
        select.addEventListener('change', function () { Units.set(quantity, select.value); });
    }
}

document.addEventListener('DOMContentLoaded', function () {
    GasLibrary.load();
    Units.load();
    setupUnitPickers();
    Units.updateLabels(document);
    Units.applyDefaults(document);

    const markStale = function () {
        if (appState.shock) byId('shockOutput').classList.add('stale');
    };
    drivenPicker = new GasPicker(byId('drivenGas'), byId('drivenMixture'), { label: 'Driven gas', onChange: markStale });
    driverPicker = new GasPicker(byId('driverGas'), byId('driverMixture'), { label: 'Driver gas', onChange: markStale });
    const defaults = { driven: { gasId: '1' }, driver: { gasId: '1' } };
    applyShockInputState(defaults);

    // Shock tube
    onSubmit('shockForm', calcShock);
    byId('shockForm').addEventListener('input', markStale);
    byId('shockForm').addEventListener('change', function (event) {
        if (event.target.name === 'mode') updateShockMode();
    });
    byId('shockForm').addEventListener('reset', function () {
        // Runs before the browser restores default values
        setTimeout(function () {
            Units.applyDefaults(byId('shockForm'));
            applyShockInputState({ mode: 'mach', gm: 'constant', ...defaults });
            clearInvalid(byId('shockForm'));
            clearShockOutputs();
            showMessages(byId('shockMessages'), []);
            try { localStorage.removeItem(STATE_STORAGE_KEY); } catch (e) { /* ignore */ }
            history.replaceState(null, '', location.pathname + location.hash);
        }, 0);
    });
    byId('copyLinkBtn').addEventListener('click', copyShockLink);

    // Redisplay results when the display units change
    Units.onChange(function () {
        if (appState.shock) displayShockResult(appState.shock);
        if (byId('pboosttank').value !== '') calcBoost();
    });

    // Timing and boost tank
    onSubmit('timingForm', calcTime);
    onReset('timingForm', TIMING_OUTPUT_IDS, 'timingMessages');
    onSubmit('boostForm', calcBoost);
    onReset('boostForm', BOOST_OUTPUT_IDS, 'boostMessages');

    // Custom gas library
    onSubmit('gasForm', submitGasForm);
    byId('gasFormEos').addEventListener('change', updateGasFormEos);
    byId('gasFormCancel').addEventListener('click', resetGasForm);
    GasLibrary.onChange(renderGasList);
    renderGasList();
    updateGasFormEos();

    // Equations dialog
    byId('showEquationsBtn').addEventListener('click', showEquationsModal);
    byId('closeEquationsBtn').addEventListener('click', closeEquationsModal);
    byId('equationsModal').addEventListener('click', function (event) {
        if (event.target === event.currentTarget) closeEquationsModal();
    });
    document.addEventListener('keydown', function (event) {
        if (event.key === 'Escape' && byId('equationsModal').classList.contains('open')) closeEquationsModal();
    });

    // Open a collapsed tool when it is the target of a navigation link
    const openTarget = function () {
        if (!/^#[a-z]+$/.test(location.hash)) return;
        const target = document.querySelector('details.tool' + location.hash);
        if (target) target.open = true;
    };
    document.querySelectorAll('.page-nav a').forEach(function (link) {
        link.addEventListener('click', function () {
            const target = document.querySelector('details' + link.getAttribute('href'));
            if (target) target.open = true;
        });
    });
    openTarget();

    byId('nasaSearch').addEventListener('input', searchNasa);
    byId('nasaSearch').addEventListener('focus', function () { loadNasaDatabase().catch(function () {}); });

    // A shared link is calculated immediately so the reader sees its results
    const restore = function () {
        if (loadShockState()) calcShock();
    };
    if (linkNeedsNasaDatabase()) loadNasaDatabase().then(restore, restore);
    else restore();
});
