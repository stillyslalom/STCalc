/**
 * x-t diagram tool: slab configuration, running the Euler solver, plots and
 * text exports. HDF5 export lives in xt-export.js.
 */

let xtSlabs = [];           // [{id, picker value, pressure (Pa), temperature (K), length (m)}]
let xtSlabCounter = 0;
let xtSlabPickers = {};     // slab id -> GasPicker
let xtVisualization = null;
let xtResults = null;
let xtRunConfig = null;     // settings of the run that produced xtResults
let xtProbePositions = [];

function showXTMessage(message, type) {
    showMessages(byId('xtMessageContainer'), message ? [{ type: type || 'info', text: message }] : []);
}

/* ------------------------------------------------------------------ */
/* Slab list                                                           */
/* ------------------------------------------------------------------ */

function addXTSlab(config) {
    const previous = xtSlabs[xtSlabs.length - 1];
    const c = config || {};
    xtSlabs.push({
        id: xtSlabCounter++,
        gas: c.gas || { gasId: '', mixture: [] },
        pressure: c.pressure || (previous ? previous.pressure : 101325),
        temperature: c.temperature || (previous ? previous.temperature : 300),
        length: c.length || 1
    });
}

/** Copy edits out of the gas pickers before the list is rebuilt. */
function syncXTSlabGases() {
    for (const slab of xtSlabs) {
        if (xtSlabPickers[slab.id]) slab.gas = xtSlabPickers[slab.id].getValue();
    }
}

function renderXTSlabs() {
    const container = byId('xtSlabList');
    Object.values(xtSlabPickers).forEach(function (picker) { picker.destroy(); });
    xtSlabPickers = {};

    if (xtSlabs.length === 0) {
        container.replaceChildren(el('p', { class: 'hint' }, 'No slabs.'));
        return;
    }

    const move = function (index, delta) {
        syncXTSlabGases();
        const slab = xtSlabs.splice(index, 1)[0];
        xtSlabs.splice(index + delta, 0, slab);
        renderXTSlabs();
    };

    container.replaceChildren(...xtSlabs.map(function (slab, index) {
        const n = index + 1;
        const select = el('select', { 'aria-label': 'Slab ' + n + ' gas' });
        const mixture = el('div', { class: 'mixture', hidden: true });
        const numberField = function (label, value, onInput) {
            const input = el('input', { type: 'number', step: 'any', inputMode: 'decimal', value: value });
            input.addEventListener('input', function () { onInput(input.value === '' ? NaN : Number(input.value)); });
            return el('label', { class: 'slab-field' }, el('span', {}, label), input);
        };

        const node = el('div', { class: 'slab-item' },
            el('div', { class: 'slab-number' }, '#' + n),
            el('div', { class: 'slab-controls' },
                el('label', { class: 'slab-field' }, el('span', {}, 'Material'), select),
                numberField('Pressure (' + Units.label('pressure') + ')', Units.inputText('pressure', slab.pressure),
                    function (v) { slab.pressure = Units.fromDisplay('pressure', v); }),
                numberField('Temperature (' + Units.label('temperature') + ')', Units.inputText('temperature', slab.temperature),
                    function (v) { slab.temperature = Units.fromDisplay('temperature', v); }),
                numberField('Length (m)', slab.length, function (v) { slab.length = v; })),
            el('div', { class: 'slab-actions' },
                el('button', { type: 'button', class: 'secondary small', disabled: index === 0, 'aria-label': 'Move slab ' + n + ' toward the driver end', onclick: function () { move(index, -1); } }, '↑'),
                el('button', { type: 'button', class: 'secondary small', disabled: index === xtSlabs.length - 1, 'aria-label': 'Move slab ' + n + ' toward the driven end', onclick: function () { move(index, 1); } }, '↓'),
                el('button', {
                    type: 'button', class: 'danger small',
                    onclick: function () {
                        syncXTSlabGases();
                        xtSlabs.splice(index, 1);
                        renderXTSlabs();
                    }
                }, 'Remove')),
            mixture);

        const picker = new GasPicker(select, mixture, { label: 'Slab ' + n + ' gas', onChange: updateInterfaceMethodInfo });
        picker.setValue(slab.gas);
        xtSlabPickers[slab.id] = picker;
        return node;
    }));
    updateInterfaceMethodInfo();
}

/** Driver and driven slabs from the shock tube calculation; other slabs are kept. */
function importShockIntoXT() {
    const shock = appState.shock;
    if (!shock) {
        showXTMessage('Run the shock tube calculation first.', 'warning');
        return false;
    }
    if (!shock.reachable) {
        showXTMessage('The shock tube calculation has no driver pressure to import.', 'warning');
        return false;
    }
    syncXTSlabGases();
    const driver = { gas: driverPicker.getValue(), pressure: shock.state4.p, temperature: shock.state4.T, length: xtSlabs[0] ? xtSlabs[0].length : 3 };
    const driven = { gas: drivenPicker.getValue(), pressure: shock.state1.p, temperature: shock.state1.T, length: xtSlabs[1] ? xtSlabs[1].length : 6 };
    const rest = xtSlabs.slice(2);
    xtSlabs = [];
    addXTSlab(driver);
    addXTSlab(driven);
    xtSlabs = xtSlabs.concat(rest);
    // Use the gas model the driver pressure was calculated with
    byId('xtGasModel').value = appState.shockInputs.variableGamma ? 'nasa' : 'constant';
    renderXTSlabs();
    return true;
}

function updateInterfaceMethodInfo() {
    const select = byId('xtInterfaceMethod');
    let hasLiquidInterface = false;
    try {
        const pinfs = xtSlabs.map(function (slab) { return xtSlabPickers[slab.id].resolve().pinf; });
        hasLiquidInterface = pinfs.some(function (p) { return p !== pinfs[0]; });
    } catch (e) { /* incomplete configuration: decided when the simulation is run */ }

    const variableGamma = byId('xtGasModel').value === 'nasa';
    const needsMixed = hasLiquidInterface || variableGamma;
    if (needsMixed) select.value = 'mixed';
    select.options[1].disabled = needsMixed;
    byId('xtInterfaceMethodLabel').title = hasLiquidInterface
        ? 'Mixed cell is required for gas–liquid interfaces.'
        : (variableGamma
            ? 'Mixed cell is required for the variable-γ gas model.'
            : 'Mixed cell is conservative. Sharp is faster but about 1% off in wave strengths.');
}

/* ------------------------------------------------------------------ */
/* Running the simulation                                              */
/* ------------------------------------------------------------------ */

function readXTConfiguration() {
    const tool = byId('xt');
    clearInvalid(tool);
    if (xtSlabs.length === 0) throw userError('Add at least one slab.');

    const variableGamma = byId('xtGasModel').value === 'nasa';
    const slabs = xtSlabs.map(function (slab, index) {
        const n = index + 1;
        const inputs = byId('xtSlabList').children[index].querySelectorAll('input[type="number"]');
        const props = xtSlabPickers[slab.id].resolve();
        if (!(slab.pressure > 0)) throw new FieldError(inputs[0], 'Slab ' + n + ': pressure must be a positive number.');
        if (!(slab.temperature > 0)) throw new FieldError(inputs[1], 'Slab ' + n + ': temperature must be a positive number.');
        if (!(slab.length > 0)) throw new FieldError(inputs[2], 'Slab ' + n + ': length must be a positive number.');
        return {
            gasId: props.name, gamma: props.gamma, mw: props.mw, pinf: props.pinf,
            thermo: variableGamma && !(props.pinf > 0) ? props.thermo : null,
            pressure: slab.pressure, temperature: slab.temperature, length: slab.length
        };
    });
    const totalLength = slabs.reduce(function (sum, slab) { return sum + slab.length; }, 0);

    const probeInput = byId('xtProbes');
    const probes = probeInput.value.split(/[,;\s]+/).filter(Boolean).map(function (text) {
        const x = Number(text);
        if (!isFinite(x) || x < 0 || x > totalLength) {
            throw new FieldError(probeInput, 'Probe "' + text + '" must be between 0 and ' + formatSig(totalLength, 4) + ' m.');
        }
        return x;
    });

    return {
        slabs: slabs,
        probes: probes,
        gasModel: variableGamma ? 'nasa' : 'constant',
        nx: Math.round(readNumber(byId('xtGridPoints'), 'Grid cells', { min: 50, max: 5000 })),
        finalTime: readNumber(byId('xtSimTime'), 'Simulated time', { positive: true }) / 1000,
        xtInterval: readNumber(byId('xtSnapshotInterval'), 'Snapshot interval', { positive: true }) / 1000,
        cfl: readNumber(byId('xtCflNumber'), 'CFL number', { min: 0.05, max: 0.9 }),
        order: Number(byId('xtOrder').value),
        integrator: byId('xtIntegrator').value,
        interfaceMethod: byId('xtInterfaceMethod').value
    };
}

function runXTSimulation() {
    let config;
    try {
        config = readXTConfiguration();
        if (config.finalTime / config.xtInterval > 5000) {
            throw new FieldError(byId('xtSnapshotInterval'), 'Too many snapshots (over 5000). Increase the interval.');
        }
    } catch (error) {
        reportError(byId('xtMessageContainer'), error);
        return;
    }

    const runBtn = byId('xtRunBtn');
    const progress = byId('xtProgressContainer');
    const fill = byId('xtProgressFill');
    runBtn.disabled = true;
    runBtn.textContent = 'Running…';
    progress.hidden = false;
    fill.style.width = '0%';
    fill.textContent = '0%';
    showXTMessage('');

    const finish = function () {
        runBtn.disabled = false;
        runBtn.textContent = 'Run';
        progress.hidden = true;
    };
    const fail = function (error) {
        console.error(error);
        showXTMessage('Simulation failed: ' + error.message, 'error');
        finish();
    };

    let solver;
    try {
        solver = new EulerSolver(config);
        solver.initialize(config.slabs);
    } catch (error) {
        fail(error);
        return;
    }

    // Run in short slices so the page stays responsive and the bar can repaint
    const slice = function () {
        let done;
        try {
            done = solver.advance(50);
        } catch (error) {
            fail(error);
            return;
        }
        const percent = Math.min(100, Math.round(solver.t / solver.finalTime * 100)) + '%';
        fill.style.width = percent;
        fill.textContent = percent;
        if (!done) {
            setTimeout(slice, 0);
            return;
        }

        try {
            xtResults = solver.getResults();
            xtRunConfig = config;
            xtProbePositions = config.probes;

            byId('xtCanvasContainer').hidden = false;
            if (!xtVisualization) xtVisualization = new XTVisualization('xtCanvas');
            xtVisualization.setProbes(xtProbePositions);
            xtVisualization.field = byId('xtField').value;
            xtVisualization.useLogScale = byId('xtLogScale').checked;
            xtVisualization.setData(xtResults);
            xtVisualization.render();
            renderXTProbes();
            displayXTStatistics(xtResults);
            finish();
        } catch (error) {
            fail(error);
        }
    };
    setTimeout(slice, 0);
}

function displayXTStatistics(results) {
    const stat = function (label, value) { return el('div', {}, el('dt', {}, label), el('dd', {}, value)); };
    byId('xtStatsContent').replaceChildren(
        stat('Time steps', results.timeSteps.toLocaleString()),
        stat('Simulated time', (results.t * 1000).toFixed(2) + ' ms'),
        stat('Wall time', (results.wallTime || 0).toFixed(2) + ' s'),
        stat('Snapshots', String(results.xtData.length)),
        stat('Tube length', formatSig(results.x.length * (results.x[1] - results.x[0]), 4) + ' m'));
}

function clearXTResults() {
    xtResults = null;
    xtRunConfig = null;
    byId('xtCanvasContainer').hidden = true;
    showXTMessage('');
}

/* ------------------------------------------------------------------ */
/* Probe traces                                                        */
/* ------------------------------------------------------------------ */

const PROBE_COLORS = ['#0066cc', '#d9480f', '#1e7e34', '#862e9c', '#b98900', '#0b7285'];

/** History of one field at one position: [{t, value}] from the stored snapshots. */
function probeSeries(results, x, field) {
    const dx = results.x[1] - results.x[0];
    const i = Math.max(0, Math.min(results.x.length - 1, Math.floor(x / dx)));
    return results.xtData.map(function (snapshot) { return { t: snapshot.t, value: snapshot[field][i] }; });
}

function renderXTProbes() {
    const hasProbes = xtProbePositions.length > 0;
    byId('xtProbeContainer').hidden = !hasProbes;
    byId('xtProbeExportBtn').hidden = !hasProbes;
    if (!hasProbes || !xtResults) return;

    const fieldKey = byId('xtField').value;
    const field = XT_FIELDS[fieldKey];
    const fit = fitCanvas(byId('xtProbeCanvas'), 900, 0.28, 200);
    const ctx = fit.ctx, width = fit.width, height = fit.height;
    const narrow = width < 560;
    // The legend sits beside the plot when there is room, otherwise above it
    const margin = narrow
        ? { left: 52, right: 12, top: 16 + 18 * xtProbePositions.length, bottom: 46 }
        : { left: 64, right: 130, top: 16, bottom: 46 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    const series = xtProbePositions.map(function (x) {
        return probeSeries(xtResults, x, fieldKey).map(function (point) {
            return { t: point.t * 1000, value: xtFieldToDisplay(fieldKey, point.value) };
        });
    });
    let lo = Infinity, hi = -Infinity;
    series.forEach(function (s) { s.forEach(function (point) { lo = Math.min(lo, point.value); hi = Math.max(hi, point.value); }); });
    if (!(hi > lo)) { hi = lo + 1; lo = lo - 1; }
    const tMax = series[0][series[0].length - 1].t;
    const yTicks = niceTicks(lo, hi, 5);
    const xTicks = niceTicks(0, tMax, narrow ? 4 : 8);
    const yMin = Math.min(lo, yTicks[0]), yMax = Math.max(hi, yTicks[yTicks.length - 1]);
    const sx = function (t) { return margin.left + t / tMax * plotWidth; };
    const sy = function (v) { return margin.top + plotHeight - (v - yMin) / (yMax - yMin) * plotHeight; };

    ctx.fillStyle = '#fff';
    ctx.fillRect(0, 0, width, height);
    ctx.font = '12px Arial';
    ctx.lineWidth = 1;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (const t of xTicks) {
        ctx.strokeStyle = '#e0e0e0';
        ctx.beginPath(); ctx.moveTo(sx(t), margin.top); ctx.lineTo(sx(t), margin.top + plotHeight); ctx.stroke();
        ctx.fillStyle = '#000';
        ctx.fillText(tickLabel(t, xTicks), sx(t), margin.top + plotHeight + 6);
    }
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const v of yTicks) {
        ctx.strokeStyle = '#e0e0e0';
        ctx.beginPath(); ctx.moveTo(margin.left, sy(v)); ctx.lineTo(margin.left + plotWidth, sy(v)); ctx.stroke();
        ctx.fillStyle = '#000';
        ctx.fillText(tickLabel(v, yTicks), margin.left - 8, sy(v));
    }
    ctx.strokeStyle = '#000';
    ctx.lineWidth = 1.5;
    ctx.strokeRect(margin.left, margin.top, plotWidth, plotHeight);

    series.forEach(function (s, n) {
        const color = PROBE_COLORS[n % PROBE_COLORS.length];
        ctx.strokeStyle = color;
        ctx.lineWidth = 2;
        ctx.beginPath();
        s.forEach(function (point, i) {
            if (i === 0) ctx.moveTo(sx(point.t), sy(point.value)); else ctx.lineTo(sx(point.t), sy(point.value));
        });
        ctx.stroke();

        // Legend entry
        const lx = narrow ? margin.left : width - margin.right + 12;
        const ly = narrow ? 14 + n * 18 : margin.top + 10 + n * 20;
        ctx.beginPath(); ctx.moveTo(lx, ly); ctx.lineTo(lx + 22, ly); ctx.stroke();
        ctx.fillStyle = '#000';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'middle';
        ctx.fillText('x = ' + formatSig(xtProbePositions[n], 4) + ' m', lx + 28, ly);
    });

    ctx.fillStyle = '#000';
    ctx.font = '14px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'alphabetic';
    ctx.fillText('Time (ms)', margin.left + plotWidth / 2, height - 10);
    const axisLabel = field.label + ' (' + xtFieldUnit(fieldKey) + ')';
    if (narrow) {
        // No room beside the tick labels: put the axis label in the top corner
        ctx.textAlign = 'right';
        ctx.fillText(axisLabel, width - margin.right, 18);
    } else {
        ctx.save();
        ctx.translate(16, margin.top + plotHeight / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText(axisLabel, 0, 0);
        ctx.restore();
    }
}

/* ------------------------------------------------------------------ */
/* Exports                                                             */
/* ------------------------------------------------------------------ */

function getTimestampString() {
    const now = new Date();
    const pad = function (n) { return String(n).padStart(2, '0'); };
    return now.getFullYear() + pad(now.getMonth() + 1) + pad(now.getDate()) + '_' +
           pad(now.getHours()) + pad(now.getMinutes()) + pad(now.getSeconds());
}

function csvCell(text) {
    const s = String(text);
    return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
}

function exportXTDataCSV() {
    const lines = ['Time (ms),Position (m),Density (kg/m³),Velocity (m/s),Pressure (kPa),Temperature (K),Gamma,MW (g/mol),p_inf (Pa),Material'];
    for (const snapshot of xtResults.xtData) {
        const t = (snapshot.t * 1000).toFixed(4);
        for (let i = 0; i < snapshot.p.length; i++) {
            lines.push([
                t, xtResults.x[i].toFixed(4), snapshot.rho[i].toPrecision(7), snapshot.u[i].toFixed(4),
                (snapshot.p[i] / 1000).toPrecision(7), snapshot.T[i].toFixed(4), snapshot.gamma[i].toFixed(6),
                snapshot.mw[i].toFixed(4), snapshot.pinf[i].toPrecision(6), csvCell(snapshot.gasId[i])
            ].join(','));
        }
    }
    downloadBlob(lines.join('\n') + '\n', 'text/csv', 'xtcalc_' + getTimestampString() + '_xtdata.csv');
}

function exportXTTracerCSV() {
    const lines = ['Interface,Time (s),Position (m),Velocity (m/s)'];
    xtResults.tracers.forEach(function (tracer, i) {
        for (const point of tracer.trajectory) lines.push([i + 1, point.t, point.x, point.v].join(','));
    });
    downloadBlob(lines.join('\n') + '\n', 'text/csv', 'xtcalc_' + getTimestampString() + '_tracers.csv');
}

function exportXTProbeCSV() {
    // Files always use the same units, whatever the page displays
    const fields = [['p', 'Pressure (kPa)', 1e-3], ['rho', 'Density (kg/m³)', 1], ['T', 'Temperature (K)', 1], ['u', 'Velocity (m/s)', 1]];
    const header = ['Time (ms)'];
    const columns = [];
    xtProbePositions.forEach(function (x) {
        fields.forEach(function (f) {
            header.push(csvCell(f[1] + ' at x = ' + x + ' m'));
            columns.push(probeSeries(xtResults, x, f[0]).map(function (point) { return point.value * f[2]; }));
        });
    });
    const lines = [header.join(',')];
    xtResults.xtData.forEach(function (snapshot, j) {
        lines.push([(snapshot.t * 1000).toFixed(5)].concat(columns.map(function (c) { return c[j].toPrecision(7); })).join(','));
    });
    downloadBlob(lines.join('\n') + '\n', 'text/csv', 'xtcalc_' + getTimestampString() + '_probes.csv');
}

function exportXTJSON() {
    // Typed arrays would serialize as objects keyed by index
    const plain = function (key, value) { return ArrayBuffer.isView(value) ? Array.from(value) : value; };
    const data = { configuration: xtRunConfig, results: xtResults };
    downloadBlob(JSON.stringify(data, plain), 'application/json', 'xtcalc_' + getTimestampString() + '.json');
}

const XT_EXPORTS = {
    png: function () { xtVisualization.exportImage('xtcalc_' + getTimestampString() + '.png'); },
    csv: exportXTDataCSV,
    tracers: exportXTTracerCSV,
    probes: exportXTProbeCSV,
    json: exportXTJSON,
    hdf5: function () { exportXTDataHDF5(); }
};

/* ------------------------------------------------------------------ */
/* Start-up                                                            */
/* ------------------------------------------------------------------ */

document.addEventListener('DOMContentLoaded', function () {
    const initialize = function () {
        if (xtSlabs.length > 0) return;
        // Start from the shock tube calculation when there is one
        if (!(appState.shock && appState.shock.reachable && importShockIntoXT())) {
            addXTSlab({ gas: { gasId: '4' }, pressure: 400000, temperature: 300, length: 3 });   // helium driver
            addXTSlab({ gas: { gasId: '1' }, pressure: 101325, temperature: 300, length: 6 });   // air driven
            renderXTSlabs();
        }
    };
    byId('xt').addEventListener('toggle', function () {
        if (byId('xt').open) initialize();
    });
    if (byId('xt').open) initialize();

    byId('xtImportBtn').addEventListener('click', function () {
        if (importShockIntoXT()) showXTMessage('');
    });
    byId('xtAddSlabBtn').addEventListener('click', function () {
        syncXTSlabGases();
        addXTSlab();
        renderXTSlabs();
    });
    byId('xtInterfaceMethod').addEventListener('change', updateInterfaceMethodInfo);
    byId('xtGasModel').addEventListener('change', updateInterfaceMethodInfo);
    byId('xtIntegrator').addEventListener('change', function () {
        byId('xtCflNumber').value = byId('xtIntegrator').value === 'SSP' ? 0.8 : 0.4;
    });
    byId('xtRunBtn').addEventListener('click', runXTSimulation);
    byId('xtClearBtn').addEventListener('click', clearXTResults);

    byId('xtField').addEventListener('change', function () {
        if (!xtVisualization) return;
        xtVisualization.setField(byId('xtField').value);
        renderXTProbes();
    });
    byId('xtLogScale').addEventListener('change', function () {
        if (xtVisualization) xtVisualization.setLogScale(byId('xtLogScale').checked);
    });

    Units.onChange(function () {
        if (xtSlabs.length > 0) {
            syncXTSlabGases();
            renderXTSlabs();
        }
        if (xtResults && !byId('xtCanvasContainer').hidden) {
            xtVisualization.render();
            renderXTProbes();
        }
    });

    onResizeEnd(function () {
        if (!xtResults || byId('xtCanvasContainer').hidden) return;
        xtVisualization.render();
        renderXTProbes();
    });

    byId('xtCanvasContainer').addEventListener('click', function (event) {
        const kind = event.target.dataset && event.target.dataset.export;
        if (!kind) return;
        if (!xtResults) {
            showXTMessage('No simulation results to export.', 'warning');
            return;
        }
        XT_EXPORTS[kind]();
    });
});
