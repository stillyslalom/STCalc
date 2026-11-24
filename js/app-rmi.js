/**
 * Shocked-interface (Richtmyer-Meshkov) calculator.
 * Uses the incident shock from the shock tube calculation (appState.shock).
 */

let interfacePicker;
let rmiGrowth = null;   // last plotted series, kept for redrawing on resize

function calculateInterface() {
    const form = byId('interfaceForm');
    const messages = byId('interfaceMessages');
    const results = byId('interfaceResults');
    clearInvalid(form);
    results.hidden = true;

    try {
        const shock = appState.shock, inputs = appState.shockInputs;
        if (!shock) throw userError('Run the shock tube calculation first.');

        const downstream = interfacePicker.resolve();
        const tDownstream = readNumber(byId('interfaceHeavyTemp'), 'Downstream temperature', { positive: true });
        const lambda = readNumber(byId('perturbWavelength'), 'Wavelength', { positive: true }) / 1000;
        const a0 = readNumber(byId('perturbAmplitude'), 'Initial amplitude', { positive: true }) / 1000;
        const distance = readNumber(byId('rmiDownstreamDist'), 'Diagnostic distance', { positive: true }) / 1000;

        // Left: shocked driven gas moving at u2. Right: downstream material at rest.
        const result = solveShockedInterface(
            { ...inputs.driven, p: shock.state2.p, T: shock.state2.T, u: shock.state2.u },
            { ...downstream, p: inputs.p1, T: tDownstream, u: 0 });

        const set = function (id, value) { byId(id).value = value; };
        set('interfaceVelocity', formatSig(result.interfaceVelocity));
        set('interfacePressure', formatSig(result.interfacePressure / 1e6));
        set('atwoodNumber', formatSig(result.atwoodNumber));
        set('waveConfigLeft', result.waveConfig.left === 'shock' ? 'Shock' : 'Rarefaction');
        set('waveConfigRight', result.waveConfig.right === 'shock' ? 'Shock' : 'Rarefaction');

        set('lightGasPressure', formatSig(result.lightGasPost.p / 1e6));
        set('lightGasTemp', formatSig(result.lightGasPost.T));
        set('lightGasDensity', formatSig(result.lightGasPost.rho));
        set('lightGasSound', formatSig(result.lightGasPost.a));
        set('heavyGasPressure', formatSig(result.heavyGasPost.p / 1e6));
        set('heavyGasTemp', formatSig(result.heavyGasPost.T));
        set('heavyGasDensity', formatSig(result.heavyGasPost.rho));
        set('heavyGasSound', formatSig(result.heavyGasPost.a));

        const notes = [];
        if (!result.converged) {
            notes.push({ type: 'warning', text: 'The Riemann solver did not fully converge; values are approximate.' });
        }
        if (!(result.interfaceVelocity > 0)) {
            throw userError('The interface does not move downstream, so it never reaches the diagnostic.');
        }

        const deltaU = calculateVelocityJump(result.interfaceVelocity, 0);
        const tArrival = distance / result.interfaceVelocity;
        const k = 2 * Math.PI / lambda;
        const params = calculateRMIParameters({ eta0: a0, lambda: lambda, A: result.atwoodNumber, deltaU: deltaU });

        byId('rmiWavenumber').textContent = formatSig(params.wavenumber, 4);
        byId('rmiKa0').textContent = formatSig(params.ka0, 4);
        byId('rmiGrowthRate').textContent = formatSig(params.initialGrowthRate, 4);
        byId('rmiDoublingTime').textContent = formatSig(params.doublingTime * 1000, 4);
        byId('rmiTimeAtLocation').textContent = formatSig(tArrival * 1000, 4);
        byId('rmiAmplitudeLinear').textContent =
            formatSig(rmiLinearAmplitude(a0, result.atwoodNumber, k, deltaU, tArrival) * 1000, 3);
        byId('rmiAmplitudeNonlinear').textContent =
            formatSig(rmiDimonteAmplitude(a0, result.atwoodNumber, k, deltaU, tArrival) * 1000, 3);

        if (params.ka0 > 1) {
            notes.push({ type: 'warning', text: 'ka₀ = ' + formatSig(params.ka0, 3) + ' is outside the small-amplitude range of the growth models.' });
        }
        if (result.atwoodNumber < 0) {
            notes.push({ type: 'info', text: 'Negative Atwood number: the perturbation reverses phase before growing.' });
        }

        results.hidden = false;
        rmiGrowth = calculateGrowthTimeSeries({
            eta0: a0, lambda: lambda, A: result.atwoodNumber, deltaU: deltaU, tMax: tArrival, numPoints: 200
        });
        plotRMIGrowth(rmiGrowth);
        showMessages(messages, notes);
    } catch (error) {
        reportError(messages, error);
    }
}

function plotRMIGrowth(growth) {
    const fit = fitCanvas(byId('rmiGrowthCanvas'), 620, 0.45, 220);
    const ctx = fit.ctx, width = fit.width, height = fit.height;
    const margin = { left: 58, right: 16, top: 16, bottom: 46 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    const times = growth.times.map(function (t) { return t * 1000; });
    const series = [
        { data: growth.linear.map(function (a) { return a * 1000; }), color: '#0066cc' },
        { data: growth.dimonte.map(function (a) { return a * 1000; }), color: '#1e7e34' }
    ];

    const all = series[0].data.concat(series[1].data);
    const yTicks = niceTicks(Math.min(0, ...all), Math.max(...all), 6);
    const xTicks = niceTicks(0, times[times.length - 1], width < 500 ? 4 : 6);
    const yMin = yTicks[0], yMax = Math.max(yTicks[yTicks.length - 1], Math.max(...all));
    const xMax = times[times.length - 1];
    const sx = function (t) { return margin.left + t / xMax * plotWidth; };
    const sy = function (a) { return margin.top + plotHeight - (a - yMin) / (yMax - yMin) * plotHeight; };

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
    for (const a of yTicks) {
        ctx.strokeStyle = '#e0e0e0';
        ctx.beginPath(); ctx.moveTo(margin.left, sy(a)); ctx.lineTo(margin.left + plotWidth, sy(a)); ctx.stroke();
        ctx.fillStyle = '#000';
        ctx.fillText(tickLabel(a, yTicks), margin.left - 8, sy(a));
    }

    ctx.strokeStyle = '#000';
    ctx.lineWidth = 1.5;
    ctx.strokeRect(margin.left, margin.top, plotWidth, plotHeight);

    for (const s of series) {
        ctx.strokeStyle = s.color;
        ctx.lineWidth = 2;
        ctx.beginPath();
        s.data.forEach(function (a, i) {
            if (i === 0) ctx.moveTo(sx(times[i]), sy(a)); else ctx.lineTo(sx(times[i]), sy(a));
        });
        ctx.stroke();
    }

    ctx.fillStyle = '#000';
    ctx.font = '14px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'alphabetic';
    ctx.fillText('Time (ms)', margin.left + plotWidth / 2, height - 10);
    ctx.save();
    ctx.translate(16, margin.top + plotHeight / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('Amplitude (mm)', 0, 0);
    ctx.restore();
}

document.addEventListener('DOMContentLoaded', function () {
    interfacePicker = new GasPicker(byId('interfaceHeavyGas'), byId('interfaceHeavyMixture'), { label: 'Downstream gas' });
    interfacePicker.setValue({ gasId: '8' });   // SF6: the classic air/SF6 experiment

    byId('interfaceForm').addEventListener('submit', function (event) {
        event.preventDefault();
        calculateInterface();
    });

    onResizeEnd(function () {
        if (rmiGrowth && !byId('interfaceResults').hidden) plotRMIGrowth(rmiGrowth);
    });

    // Results describe the previous incident shock once the shock tube inputs change
    appState.shockListeners.push(function () {
        if (!byId('interfaceResults').hidden) {
            byId('interfaceResults').hidden = true;
            rmiGrowth = null;
            showMessages(byId('interfaceMessages'), [{ type: 'info', text: 'Shock tube inputs changed — recalculate.' }]);
        }
    });
});
