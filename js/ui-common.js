/**
 * Shared UI building blocks: DOM helpers, inline messages, number parsing,
 * the custom-gas library and the gas picker used by every calculator.
 */

/* ------------------------------------------------------------------ */
/* DOM helpers                                                         */
/* ------------------------------------------------------------------ */

function byId(id) {
    return document.getElementById(id);
}

/**
 * Create an element. Text children are inserted as text nodes, so user
 * supplied strings (gas names) can never be interpreted as markup.
 */
function el(tag, props, ...children) {
    const node = document.createElement(tag);
    for (const key in props || {}) {
        const value = props[key];
        if (value === undefined || value === null || value === false) continue;
        if (key === 'class') node.className = value;
        else if (key === 'dataset') Object.assign(node.dataset, value);
        else if (key.slice(0, 2) === 'on') node.addEventListener(key.slice(2), value);
        else if (key in node) node[key] = value;
        else node.setAttribute(key, value);
    }
    for (const child of children) {
        if (child === null || child === undefined || child === false) continue;
        node.append(child);
    }
    return node;
}

/** Significant-figure formatting for output fields; '—' for missing values. */
function formatSig(value, digits = 6) {
    if (value === null || value === undefined || !isFinite(value)) return '—';
    return Number(value).toPrecision(digits);
}

/**
 * Show a list of messages in a container, replacing whatever was there.
 * @param {HTMLElement} container
 * @param {Array<{type:'error'|'warning'|'info'|'success', text:string}>} items
 */
function showMessages(container, items) {
    container.replaceChildren(...(items || []).map(function (item) {
        return el('div', { class: 'message ' + item.type }, item.text);
    }));
}

class FieldError extends Error {
    constructor(input, message) {
        super(message);
        this.input = input;
    }
}

/**
 * Read a number from an input element.
 * Unlike parseFloat this rejects partially numeric text such as "2,5" or "1.5x".
 * @throws {FieldError}
 */
function readNumber(input, label, options) {
    const opts = options || {};
    const text = String(input.value).trim();
    if (text === '') {
        if (opts.optional) return null;
        throw new FieldError(input, label + ' is required' + (input.validity && input.validity.badInput ? ' and must be a number (use a decimal point).' : '.'));
    }
    const value = Number(text);
    if (!isFinite(value)) throw new FieldError(input, label + ' must be a number (use a decimal point).');
    if (opts.positive && value <= 0) throw new FieldError(input, label + ' must be greater than zero.');
    if (opts.min !== undefined && value < opts.min) throw new FieldError(input, label + ' must be at least ' + opts.min + '.');
    if (opts.max !== undefined && value > opts.max) throw new FieldError(input, label + ' must be at most ' + opts.max + '.');
    return value;
}

function clearInvalid(root) {
    root.querySelectorAll('.invalid').forEach(function (node) {
        node.classList.remove('invalid');
        node.removeAttribute('aria-invalid');
    });
}

function markInvalid(input) {
    if (!input || !input.classList) return;
    input.classList.add('invalid');
    input.setAttribute('aria-invalid', 'true');
    if (typeof input.focus === 'function') input.focus();
}

/** Report a thrown error in a message container; marks the offending field. */
function reportError(container, error) {
    if (error instanceof FieldError) markInvalid(error.input);
    else if (!(error instanceof ShockInputError) && !error.userFacing) console.error(error);
    showMessages(container, [{ type: 'error', text: error.message }]);
}

/** An error whose message is meant to be shown to the user as is. */
function userError(message, input) {
    const error = input ? new FieldError(input, message) : new Error(message);
    error.userFacing = true;
    return error;
}

/** "Nice" axis ticks (multiples of 1, 2 or 5) covering [min, max]. */
function niceTicks(min, max, target) {
    if (!(max > min)) return [min];
    const rough = (max - min) / Math.max(1, target || 5);
    const power = Math.pow(10, Math.floor(Math.log10(rough)));
    const fraction = rough / power;
    const step = (fraction < 1.5 ? 1 : fraction < 3.5 ? 2 : fraction < 7.5 ? 5 : 10) * power;
    const ticks = [];
    for (let k = Math.ceil(min / step - 1e-9); k * step <= max + step * 1e-9; k++) {
        ticks.push(Math.abs(k * step) < step * 1e-9 ? 0 : k * step);
    }
    return ticks;
}

/** Label for a tick value with only as many decimals as the spacing needs. */
function tickLabel(value, ticks) {
    const step = ticks.length > 1 ? Math.abs(ticks[1] - ticks[0]) : Math.abs(value) || 1;
    const decimals = Math.max(0, Math.min(6, -Math.floor(Math.log10(step) + 1e-9)));
    return value.toFixed(decimals);
}

function downloadBlob(content, mimeType, filename) {
    const url = URL.createObjectURL(content instanceof Blob ? content : new Blob([content], { type: mimeType }));
    const link = el('a', { href: url, download: filename });
    link.click();
    URL.revokeObjectURL(url);
}

/* ------------------------------------------------------------------ */
/* Custom gas library (persisted in localStorage)                      */
/* ------------------------------------------------------------------ */

const GasLibrary = {
    STORAGE_KEY: 'wistl_custom_gases',
    custom: {},
    listeners: [],
    storageWorks: true,

    load() {
        this.custom = {};
        try {
            const stored = JSON.parse(localStorage.getItem(this.STORAGE_KEY) || '[]');
            for (const gas of stored) {
                if (!gas || typeof gas.id !== 'string' || validateGasDefinition(gas)) continue;
                this.custom[gas.id] = { id: gas.id, ...normalizeGas(gas) };
            }
        } catch (e) {
            this.storageWorks = false;
            console.warn('Could not read custom gases from browser storage.', e);
        }
    },

    save() {
        try {
            localStorage.setItem(this.STORAGE_KEY, JSON.stringify(Object.values(this.custom)));
            this.storageWorks = true;
        } catch (e) {
            this.storageWorks = false;
        }
        this.listeners.forEach(function (listener) { listener(); });
    },

    /** All gases keyed by id: built-in first, then custom. */
    all() {
        const gases = {};
        for (const id in GAS_DATABASE) gases[id] = normalizeGas(GAS_DATABASE[id]);
        for (const id in this.custom) gases[id] = this.custom[id];
        return gases;
    },

    get(id) {
        return this.all()[id] || null;
    },

    isCustom(id) {
        return Object.prototype.hasOwnProperty.call(this.custom, id);
    },

    findByName(name, exceptId) {
        const wanted = name.trim().toLowerCase();
        const gases = this.all();
        for (const id in gases) {
            if (id !== exceptId && gases[id].name.toLowerCase() === wanted) return id;
        }
        return null;
    },

    nextId() {
        let n = 0;
        while (this.custom['custom_' + n]) n++;
        return 'custom_' + n;
    },

    /** Add or update a gas. Returns its id. Throws a user-facing Error. */
    put(def, id) {
        const gas = normalizeGas({ ...def, name: String(def.name || '').trim() });
        const problem = validateGasDefinition(gas);
        if (problem) throw userError(problem);
        if (this.findByName(gas.name, id)) throw userError('A gas named "' + gas.name + '" already exists.');
        const gasId = id || this.nextId();
        this.custom[gasId] = { id: gasId, ...gas };
        this.save();
        return gasId;
    },

    remove(id) {
        delete this.custom[id];
        this.save();
    },

    /**
     * Id of a gas with exactly these properties, adding it if necessary.
     * Used when a shared link refers to a custom gas this browser lacks.
     */
    ensure(def) {
        const gas = normalizeGas(def);
        const gases = this.all();
        for (const id in gases) {
            const g = gases[id];
            if (g.name === gas.name && g.gamma === gas.gamma && g.mw === gas.mw && g.pinf === gas.pinf) return id;
        }
        let name = gas.name;
        for (let n = 2; this.findByName(name); n++) name = gas.name + ' (' + n + ')';
        return this.put({ ...gas, name: name });
    },

    onChange(listener) {
        this.listeners.push(listener);
    }
};

/* ------------------------------------------------------------------ */
/* Gas picker: one <select> plus an optional mole-fraction mixture     */
/* ------------------------------------------------------------------ */

const MIXTURE_VALUE = 'mixture';

function fillGasSelect(select, selected, allowMixture) {
    const gases = GasLibrary.all();
    const groups = { 'Gases': [], 'Liquids': [], 'Custom': [] };
    for (const id in gases) {
        const group = GasLibrary.isCustom(id) ? 'Custom' : (isStiffened(gases[id]) ? 'Liquids' : 'Gases');
        groups[group].push([id, gases[id].name]);
    }
    const children = [el('option', { value: '' }, 'Select…')];
    for (const label in groups) {
        if (groups[label].length === 0) continue;
        groups[label].sort(function (a, b) { return a[1].localeCompare(b[1]); });
        children.push(el('optgroup', { label: label }, ...groups[label].map(function (entry) {
            return el('option', { value: entry[0] }, entry[1]);
        })));
    }
    if (allowMixture) children.push(el('option', { value: MIXTURE_VALUE }, 'Gas mixture…'));
    select.replaceChildren(...children);
    select.value = selected || '';
    if (select.value !== (selected || '')) select.value = '';   // selected gas no longer exists
}

class GasPicker {
    /**
     * @param {HTMLSelectElement} select
     * @param {HTMLElement} mixtureContainer - shown when "Gas mixture…" is chosen
     * @param {Object} [options] - {label, onChange}
     */
    constructor(select, mixtureContainer, options) {
        this.select = select;
        this.container = mixtureContainer;
        this.options = options || {};
        this.rows = [];
        this.refresh = this.refresh.bind(this);

        fillGasSelect(select, '', true);
        select.addEventListener('change', () => {
            if (select.value === MIXTURE_VALUE && this.rows.length === 0) {
                this.rows = [{ gasId: '', fraction: '' }, { gasId: '', fraction: '' }];
            }
            this.render();
            this.changed();
        });
        GasPicker.instances.add(this);
        this.render();
    }

    /** Stop tracking this picker (its elements are about to be discarded). */
    destroy() {
        GasPicker.instances.delete(this);
    }

    changed() {
        if (this.options.onChange) this.options.onChange();
    }

    /** Rebuild option lists after the gas library changed. */
    refresh() {
        fillGasSelect(this.select, this.select.value, true);
        this.render();
    }

    getValue() {
        return {
            gasId: this.select.value,
            mixture: this.rows.map(function (row) { return { gasId: row.gasId, fraction: row.fraction }; })
        };
    }

    setValue(value) {
        const v = value || {};
        this.rows = (v.mixture || []).map(function (row) {
            return { gasId: row.gasId || '', fraction: row.fraction === undefined ? '' : String(row.fraction) };
        });
        // 'custom' was the mixture value in states saved by earlier versions
        const gasId = v.gasId === 'custom' ? MIXTURE_VALUE : (v.gasId || '');
        fillGasSelect(this.select, gasId, true);
        this.render();
    }

    total() {
        return this.rows.reduce(function (sum, row) { return sum + (parseFloat(row.fraction) || 0); }, 0);
    }

    render() {
        const isMixture = this.select.value === MIXTURE_VALUE;
        this.container.hidden = !isMixture;
        if (!isMixture) {
            this.container.replaceChildren();
            return;
        }

        const label = this.options.label || 'Gas';
        const totalNode = el('span', { class: 'mixture-total' });
        const updateTotal = () => {
            const total = this.total();
            const ok = Math.abs(total - 1) < 0.001;
            totalNode.className = 'mixture-total ' + (ok ? 'ok' : 'bad');
            totalNode.textContent = 'Total: ' + (Math.round(total * 1e4) / 100) + '%' + (ok ? ' ✓' : ' (must be 100%)');
        };

        const rowNodes = this.rows.map((row, index) => {
            const gasSelect = el('select', { 'aria-label': label + ' mixture component ' + (index + 1) });
            fillGasSelect(gasSelect, row.gasId, false);
            row.gasId = gasSelect.value;
            gasSelect.addEventListener('change', () => {
                row.gasId = gasSelect.value;
                this.changed();
            });
            const fraction = el('input', {
                type: 'number', step: 'any', min: '0', max: '1', inputMode: 'decimal', value: row.fraction,
                placeholder: '0–1', 'aria-label': label + ' mole fraction of component ' + (index + 1)
            });
            fraction.addEventListener('input', () => {
                row.fraction = fraction.value;
                updateTotal();
                this.changed();
            });
            const remove = el('button', {
                type: 'button', class: 'danger small', disabled: this.rows.length <= 1,
                onclick: () => {
                    this.rows.splice(index, 1);
                    this.render();
                    this.changed();
                }
            }, 'Remove');
            return el('div', { class: 'mixture-row' }, gasSelect, el('span', {}, 'mole fraction'), fraction, remove);
        });

        const add = el('button', {
            type: 'button', class: 'secondary small',
            onclick: () => {
                this.rows.push({ gasId: '', fraction: '' });
                this.render();
                this.changed();
            }
        }, 'Add component');

        this.container.replaceChildren(...rowNodes, el('div', {}, add, totalNode));
        updateTotal();
    }

    /**
     * Thermodynamic properties of the current selection.
     * @returns {{gamma:number, mw:number, pinf:number, name:string}}
     * @throws {FieldError} when nothing valid is selected
     */
    resolve() {
        const label = this.options.label || 'Gas';
        if (this.select.value === MIXTURE_VALUE) {
            try {
                return mixtureProperties(this.rows, GasLibrary.all());
            } catch (e) {
                throw new FieldError(this.select, label + ': ' + e.message);
            }
        }
        const gas = GasLibrary.get(this.select.value);
        if (!gas) throw new FieldError(this.select, 'Select a ' + label.toLowerCase() + '.');
        return { gamma: gas.gamma, mw: gas.mw, pinf: gas.pinf, name: gas.name };
    }
}

GasPicker.instances = new Set();
GasLibrary.onChange(function () {
    GasPicker.instances.forEach(function (picker) { picker.refresh(); });
});

/* ------------------------------------------------------------------ */
/* Compact text form of a gas selection, for shareable links           */
/* ------------------------------------------------------------------ */

function encodeGasId(id) {
    if (!GasLibrary.isCustom(id)) return encodeURIComponent(id);
    const g = GasLibrary.custom[id];
    return ['c', encodeURIComponent(g.name), g.gamma, g.mw, g.pinf].join(':');
}

function decodeGasId(text) {
    const parts = text.split(':');
    if (parts[0] === 'c' && parts.length === 5) {
        return GasLibrary.ensure({
            name: decodeURIComponent(parts[1]),
            gamma: Number(parts[2]), mw: Number(parts[3]), pinf: Number(parts[4])
        });
    }
    return decodeURIComponent(text);
}

/** e.g. "1", "c:Methane:1.31:16.04:0" or "mix,4@0.9,2@0.1" */
function encodeGasValue(value) {
    if (value.gasId !== MIXTURE_VALUE) return encodeGasId(value.gasId);
    return ['mix'].concat(value.mixture.filter(function (row) { return row.gasId; }).map(function (row) {
        return encodeGasId(row.gasId) + '@' + (parseFloat(row.fraction) || 0);
    })).join(',');
}

function decodeGasValue(text) {
    if (text.slice(0, 3) !== 'mix') return { gasId: decodeGasId(text), mixture: [] };
    return {
        gasId: MIXTURE_VALUE,
        mixture: text.split(',').slice(1).map(function (part) {
            const at = part.lastIndexOf('@');
            return { gasId: decodeGasId(part.slice(0, at)), fraction: part.slice(at + 1) };
        })
    };
}
