/**
 * X-T diagram rendering: a heatmap of one flow field over position and time,
 * with material interface paths and optional probe locations.
 */

// matplotlib "inferno" colormap as consecutive RGB triples
const INFERNO_LUT = new Uint8Array([
    0,0,4,1,0,5,1,1,6,1,1,8,2,1,10,2,2,12,2,2,14,3,2,16,
    4,3,18,4,3,20,5,4,23,6,4,25,7,5,27,8,5,29,9,6,31,10,7,34,
    11,7,36,12,8,38,13,8,41,14,9,43,16,9,45,17,10,48,18,10,50,20,11,52,
    21,11,55,22,11,57,24,12,60,25,12,62,27,12,65,28,12,67,30,12,69,31,12,72,
    33,12,74,35,12,76,36,12,79,38,12,81,40,11,83,41,11,85,43,11,87,45,11,89,
    47,10,91,49,10,92,50,10,94,52,10,95,54,9,97,56,9,98,57,9,99,59,9,100,
    61,9,101,62,9,102,64,10,103,66,10,104,68,10,104,69,10,105,71,11,106,73,11,106,
    74,12,107,76,12,107,77,13,108,79,13,108,81,14,108,82,14,109,84,15,109,85,15,109,
    87,16,110,89,16,110,90,17,110,92,18,110,93,18,110,95,19,110,96,19,110,98,20,110,
    100,21,110,101,21,110,103,22,110,104,22,110,106,23,110,108,24,110,109,24,110,111,25,110,
    112,25,110,114,26,110,115,26,110,117,27,110,119,28,109,120,28,109,122,29,109,123,29,109,
    125,30,109,126,30,109,128,31,108,129,31,108,131,32,108,133,32,107,134,33,107,136,33,107,
    137,34,106,139,34,106,140,35,105,142,35,105,144,36,105,145,36,104,147,37,104,148,37,103,
    150,38,103,151,38,102,153,39,102,154,39,101,156,40,100,157,40,100,159,41,99,160,41,99,
    162,42,98,163,42,97,165,43,97,166,43,96,168,44,95,169,44,95,171,45,94,172,46,93,
    174,46,92,175,47,92,177,47,91,178,48,90,180,48,89,181,49,88,183,49,88,184,50,87,
    186,50,86,187,51,85,188,51,84,190,52,83,191,52,82,193,53,81,194,53,80,196,54,80,
    197,54,79,199,55,78,200,55,77,201,56,76,203,56,75,204,57,74,206,57,73,207,58,72,
    208,58,71,210,59,70,211,59,69,213,60,68,214,60,67,215,61,66,217,61,65,218,62,64,
    219,62,62,221,63,61,222,63,60,223,64,59,225,64,58,226,65,57,227,65,56,228,66,55,
    230,66,54,231,67,53,232,67,52,234,68,51,235,68,50,236,69,49,237,69,48,239,70,47,
    240,70,45,241,71,44,243,71,43,244,72,42,245,72,41,246,73,40,248,73,39,249,74,38,
    250,74,37,251,75,36,253,75,35,254,76,34,255,76,33,255,77,33,255,78,33,255,79,33,
    255,80,33,255,81,34,255,82,35,255,83,35,255,84,36,255,85,37,255,86,38,255,87,39,
    255,88,40,255,89,41,255,90,42,255,91,43,255,93,45,255,94,46,255,95,47,255,96,48,
    255,97,50,255,98,51,255,99,53,255,100,54,255,101,56,255,102,57,255,103,59,255,104,61,
    255,105,62,255,106,64,255,107,66,255,108,68,255,109,70,255,110,71,255,112,73,255,113,75,
    255,114,77,255,115,79,255,116,81,255,117,83,255,118,85,255,119,87,255,120,89,255,122,91,
    255,123,94,255,124,96,255,125,98,255,126,100,255,127,103,255,128,105,255,129,107,255,130,109,
    255,131,112,255,132,114,255,133,116,255,134,119,255,135,121,255,136,123,255,137,126,255,138,128,
    255,140,131,255,141,133,255,142,136,255,143,138,255,144,141,255,145,143,255,146,146,255,147,148,
    255,148,151,255,149,153,255,150,156,255,151,159,255,152,161,255,153,164,255,154,167,255,155,169,
    255,156,172,255,157,175,255,158,178,255,159,180,255,160,183,255,161,186,255,162,189,255,163,192,
    255,164,194,255,165,197,255,166,200,255,167,203,255,168,206,255,169,209,255,170,212,255,171,215,
    255,172,218,255,173,220,255,174,223,255,175,226,255,176,229,255,177,232,255,178,235,255,179,238,
    255,180,241,255,181,244,255,182,247,255,182,250,255,183,253,255,184,255
]);

const INFERNO_LAST = INFERNO_LUT.length / 3 - 1;

const XT_FIELDS = {
    p:   { label: 'Pressure', unit: 'kPa', scale: 1e-3 },
    rho: { label: 'Density', unit: 'kg/m³', scale: 1 },
    T:   { label: 'Temperature', unit: 'K', scale: 1 },
    u:   { label: 'Velocity', unit: 'm/s', scale: 1 }
};

class XTVisualization {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas.getContext('2d');
        // The finished plot is cached so that tooltips redraw cheaply
        this.cacheCanvas = document.createElement('canvas');

        this.results = null;
        this.field = 'p';
        this.useLogScale = false;
        this.probes = [];
        this.vMin = 0;
        this.vMax = 1;

        this.mouseX = -1;
        this.mouseY = -1;
        this.tooltipPending = false;
        this.layout();
        this.setupEventListeners();
    }

    /** Fit the canvas to its container; all drawing uses CSS-pixel coordinates. */
    layout() {
        const fit = fitCanvas(this.canvas, 1000, 0.6, 300);
        this.ctx = fit.ctx;
        this.width = fit.width;
        this.height = fit.height;
        this.dpr = fit.dpr;
        this.narrow = this.width < 560;
        this.margin = this.narrow
            ? { top: 30, right: 62, bottom: 44, left: 40 }
            : { top: 36, right: 84, bottom: 52, left: 58 };
        this.plotWidth = this.width - this.margin.left - this.margin.right;
        this.plotHeight = this.height - this.margin.top - this.margin.bottom;
        this.cacheCanvas.width = this.canvas.width;
        this.cacheCanvas.height = this.canvas.height;
    }

    setupEventListeners() {
        const move = (clientX, clientY) => {
            const rect = this.canvas.getBoundingClientRect();
            this.mouseX = (clientX - rect.left) * this.width / rect.width;
            this.mouseY = (clientY - rect.top) * this.height / rect.height;
            if (!this.tooltipPending) {
                this.tooltipPending = true;
                requestAnimationFrame(() => {
                    this.tooltipPending = false;
                    this.drawTooltip();
                });
            }
        };
        this.canvas.addEventListener('mousemove', e => move(e.clientX, e.clientY));
        this.canvas.addEventListener('touchstart', e => move(e.touches[0].clientX, e.touches[0].clientY), { passive: true });
        this.canvas.addEventListener('touchmove', e => move(e.touches[0].clientX, e.touches[0].clientY), { passive: true });
        this.canvas.addEventListener('mouseleave', () => {
            this.mouseX = -1;
            this.mouseY = -1;
            this.drawFromCache();
        });
    }

    /** Results are only read, never modified. */
    setData(results) {
        this.results = results;
        const x = results.x;
        this.dx = x.length > 1 ? x[1] - x[0] : 1;
        this.xMin = 0;
        this.xMax = x[x.length - 1] + 0.5 * this.dx;
        this.tMin = 0;
        this.tMax = results.xtData[results.xtData.length - 1].t;
        this.updateRange();
    }

    setField(field) {
        this.field = field;
        if (this.results) {
            this.updateRange();
            this.render();
        }
    }

    setLogScale(enabled) {
        this.useLogScale = enabled;
        this.render();
    }

    setProbes(positions) {
        this.probes = positions || [];
    }

    /** Logarithmic colour needs strictly positive data. */
    logScaleActive() {
        return this.useLogScale && this.vMin > 0 && this.vMax > this.vMin;
    }

    updateRange() {
        let min = Infinity, max = -Infinity;
        for (const snapshot of this.results.xtData) {
            const values = snapshot[this.field];
            for (let i = 0; i < values.length; i++) {
                if (values[i] < min) min = values[i];
                if (values[i] > max) max = values[i];
            }
        }
        this.vMin = min;
        this.vMax = max;
    }

    normalize(value) {
        if (!(this.vMax > this.vMin)) return 0.5;
        if (this.logScaleActive()) {
            return Math.log(value / this.vMin) / Math.log(this.vMax / this.vMin);
        }
        return (value - this.vMin) / (this.vMax - this.vMin);
    }

    xToPixel(x) {
        return this.margin.left + (x - this.xMin) / (this.xMax - this.xMin) * this.plotWidth;
    }

    tToPixel(t) {
        // Time increases upward
        return this.margin.top + this.plotHeight - (t - this.tMin) / (this.tMax - this.tMin) * this.plotHeight;
    }

    /** Snapshot index and cell index for a point in data coordinates. */
    locate(x, t) {
        const data = this.results.xtData;
        const nx = this.results.x.length;
        const j = Math.max(0, Math.min(data.length - 1, Math.round((t - this.tMin) / (this.tMax - this.tMin) * (data.length - 1))));
        const i = Math.max(0, Math.min(nx - 1, Math.floor(x / this.dx)));
        return { i, j };
    }

    render() {
        this.layout();
        const visible = this.ctx;
        this.ctx = this.cacheCanvas.getContext('2d');
        this.ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
        this.ctx.fillStyle = '#ffffff';
        this.ctx.fillRect(0, 0, this.width, this.height);
        if (this.results) {
            this.drawHeatmap();
            this.drawTracers();
            this.drawProbes();
            this.drawAxes();
            this.drawColorbar();
        }
        this.ctx = visible;
        this.drawFromCache();
    }

    drawFromCache() {
        this.ctx.save();
        this.ctx.setTransform(1, 0, 0, 1, 0, 0);   // the cache is in device pixels
        this.ctx.drawImage(this.cacheCanvas, 0, 0);
        this.ctx.restore();
    }

    drawHeatmap() {
        const data = this.results.xtData;
        const nx = this.results.x.length;
        // Image data is addressed in device pixels
        const x0 = Math.round(this.margin.left * this.dpr);
        const y0 = Math.round(this.margin.top * this.dpr);
        const w = Math.round((this.margin.left + this.plotWidth) * this.dpr) - x0;
        const h = Math.round((this.margin.top + this.plotHeight) * this.dpr) - y0;
        const image = this.ctx.createImageData(w, h);
        const pixels = image.data;

        // Nearest cell for every pixel column
        const column = new Int32Array(w);
        for (let px = 0; px < w; px++) {
            column[px] = Math.min(nx - 1, Math.floor((px + 0.5) / w * nx));
        }

        for (let py = 0; py < h; py++) {
            const fraction = 1 - (py + 0.5) / h;
            const values = data[Math.min(data.length - 1, Math.round(fraction * (data.length - 1)))][this.field];
            for (let px = 0; px < w; px++) {
                const norm = Math.max(0, Math.min(1, this.normalize(values[column[px]])));
                const c = Math.round(norm * INFERNO_LAST) * 3;
                const o = (py * w + px) * 4;
                pixels[o] = INFERNO_LUT[c];
                pixels[o + 1] = INFERNO_LUT[c + 1];
                pixels[o + 2] = INFERNO_LUT[c + 2];
                pixels[o + 3] = 255;
            }
        }
        this.ctx.putImageData(image, x0, y0);
    }

    drawTracers() {
        this.ctx.strokeStyle = '#ffffff';
        this.ctx.lineWidth = 2;
        this.ctx.setLineDash([]);
        for (const tracer of this.results.tracers || []) {
            this.ctx.beginPath();
            tracer.trajectory.forEach((point, i) => {
                const x = this.xToPixel(point.x), y = this.tToPixel(point.t);
                if (i === 0) this.ctx.moveTo(x, y); else this.ctx.lineTo(x, y);
            });
            this.ctx.stroke();
        }
    }

    drawProbes() {
        this.ctx.strokeStyle = '#35c9ff';
        this.ctx.lineWidth = 1.5;
        this.ctx.setLineDash([6, 4]);
        for (const x of this.probes) {
            this.ctx.beginPath();
            this.ctx.moveTo(this.xToPixel(x), this.margin.top);
            this.ctx.lineTo(this.xToPixel(x), this.margin.top + this.plotHeight);
            this.ctx.stroke();
        }
        this.ctx.setLineDash([]);
    }

    drawAxes() {
        const ctx = this.ctx;
        const field = XT_FIELDS[this.field];
        const bottom = this.margin.top + this.plotHeight;
        ctx.strokeStyle = '#000000';
        ctx.fillStyle = '#000000';
        ctx.lineWidth = 1.5;
        ctx.strokeRect(this.margin.left, this.margin.top, this.plotWidth, this.plotHeight);

        ctx.font = '12px Arial';
        ctx.textAlign = 'center';
        const xTicks = niceTicks(this.xMin, this.xMax, this.narrow ? 5 : 8);
        for (const x of xTicks) {
            const px = this.xToPixel(x);
            ctx.beginPath(); ctx.moveTo(px, bottom); ctx.lineTo(px, bottom + 5); ctx.stroke();
            ctx.fillText(tickLabel(x, xTicks), px, bottom + 20);
        }
        ctx.textAlign = 'right';
        const tTicks = niceTicks(this.tMin * 1000, this.tMax * 1000, 6);
        for (const t of tTicks) {
            const py = this.tToPixel(t / 1000);
            ctx.beginPath(); ctx.moveTo(this.margin.left - 5, py); ctx.lineTo(this.margin.left, py); ctx.stroke();
            ctx.fillText(tickLabel(t, tTicks), this.margin.left - 8, py + 4);
        }

        ctx.font = '13px Arial';
        ctx.textAlign = 'center';
        ctx.fillText('Position (m)', this.margin.left + this.plotWidth / 2, this.height - 8);
        ctx.save();
        ctx.translate(12, this.margin.top + this.plotHeight / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText('Time (ms)', 0, 0);
        ctx.restore();

        // The colour scale is labelled above the colour bar
        ctx.font = 'bold 13px Arial';
        ctx.textAlign = 'right';
        ctx.fillText(field.label + ' (' + field.unit + ')', this.width - 6, this.margin.top - 12);
    }

    drawColorbar() {
        const ctx = this.ctx;
        const field = XT_FIELDS[this.field];
        const barWidth = 14;
        const barX = this.width - this.margin.right + 10;
        const barY = this.margin.top;
        const barHeight = this.plotHeight;

        for (let py = 0; py < barHeight; py++) {
            const c = Math.round((1 - py / (barHeight - 1)) * INFERNO_LAST) * 3;
            ctx.fillStyle = `rgb(${INFERNO_LUT[c]},${INFERNO_LUT[c + 1]},${INFERNO_LUT[c + 2]})`;
            ctx.fillRect(barX, barY + py, barWidth, 1);
        }
        ctx.strokeStyle = '#000000';
        ctx.fillStyle = '#000000';
        ctx.lineWidth = 1;
        ctx.strokeRect(barX, barY, barWidth, barHeight);

        ctx.font = '11px Arial';
        ctx.textAlign = 'left';
        const lo = this.vMin * field.scale, hi = this.vMax * field.scale;
        let ticks;
        if (this.logScaleActive()) {
            ticks = [];
            for (let e = Math.ceil(Math.log10(lo) - 1e-9); e <= Math.floor(Math.log10(hi) + 1e-9); e++) ticks.push(Math.pow(10, e));
            if (ticks.length < 2) ticks = [lo, Math.sqrt(lo * hi), hi];
        } else {
            ticks = niceTicks(lo, hi, this.narrow ? 4 : 6);
        }
        for (const value of ticks) {
            const norm = this.normalize(value / field.scale);
            if (norm < -1e-9 || norm > 1 + 1e-9) continue;
            const y = barY + barHeight - norm * barHeight;
            ctx.beginPath(); ctx.moveTo(barX + barWidth, y); ctx.lineTo(barX + barWidth + 4, y); ctx.stroke();
            const label = this.logScaleActive() ? Number(value.toPrecision(3)).toString() : tickLabel(value, ticks);
            ctx.fillText(label, barX + barWidth + 7, y + 4);
        }

    }

    drawTooltip() {
        if (!this.results) return;
        this.drawFromCache();
        const inside = this.mouseX >= this.margin.left && this.mouseX <= this.margin.left + this.plotWidth &&
                       this.mouseY >= this.margin.top && this.mouseY <= this.margin.top + this.plotHeight;
        if (!inside) return;

        const ctx = this.ctx;
        const x = this.xMin + (this.mouseX - this.margin.left) / this.plotWidth * (this.xMax - this.xMin);
        const t = this.tMax - (this.mouseY - this.margin.top) / this.plotHeight * (this.tMax - this.tMin);
        const { i, j } = this.locate(x, t);
        const snapshot = this.results.xtData[j];
        const lines = [
            `x ${x.toFixed(3)} m, t ${(t * 1000).toFixed(3)} ms`,
            `p ${(snapshot.p[i] / 1000).toPrecision(5)} kPa`,
            `T ${snapshot.T[i].toPrecision(5)} K`,
            `ρ ${snapshot.rho[i].toPrecision(5)} kg/m³`,
            `u ${snapshot.u[i].toFixed(2)} m/s`
        ];

        // Crosshair
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.6)';
        ctx.lineWidth = 1;
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(this.margin.left, this.mouseY);
        ctx.lineTo(this.margin.left + this.plotWidth, this.mouseY);
        ctx.moveTo(this.mouseX, this.margin.top);
        ctx.lineTo(this.mouseX, this.margin.top + this.plotHeight);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.font = '12px Arial';
        const boxWidth = Math.max(...lines.map(line => ctx.measureText(line).width)) + 12;
        const boxHeight = lines.length * 16 + 8;
        let boxX = this.mouseX + 12;
        let boxY = this.mouseY - boxHeight - 8;
        if (boxX + boxWidth > this.width - 4) boxX = Math.max(4, this.mouseX - boxWidth - 12);
        if (boxY < this.margin.top) boxY = this.mouseY + 12;

        ctx.fillStyle = 'rgba(255, 255, 255, 0.93)';
        ctx.fillRect(boxX, boxY, boxWidth, boxHeight);
        ctx.strokeStyle = '#000000';
        ctx.strokeRect(boxX, boxY, boxWidth, boxHeight);
        ctx.fillStyle = '#000000';
        ctx.textAlign = 'left';
        lines.forEach((line, n) => ctx.fillText(line, boxX + 6, boxY + 16 * (n + 1)));
    }

    exportImage(filename = 'xt-diagram.png') {
        this.drawFromCache();
        const link = document.createElement('a');
        link.download = filename;
        link.href = this.canvas.toDataURL('image/png');
        link.click();
    }
}
