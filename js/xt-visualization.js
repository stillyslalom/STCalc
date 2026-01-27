/**
 * X-T Diagram Visualization with Plasma Colormap
 * Renders pressure heatmap and interface trajectories
 */

class XTVisualization {
    constructor(canvasId, config = {}) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas.getContext('2d');
        
        // Configuration
        this.width = config.width || 1000;
        this.height = config.height || 600;
        this.margin = config.margin || { top: 40, right: 100, bottom: 60, left: 70 };
        
        // Set canvas size
        this.canvas.width = this.width;
        this.canvas.height = this.height;
        
        // Plot area dimensions
        this.plotWidth = this.width - this.margin.left - this.margin.right;
        this.plotHeight = this.height - this.margin.top - this.margin.bottom;
        
        // Data
        this.xtData = null;
        this.tracers = null;
        this.xMin = 0;
        this.xMax = 1;
        this.tMin = 0;
        this.tMax = 1;
        this.pMin = 0;
        this.pMax = 1;
        
        // Colormap settings
        this.useLogScale = false;
        this.colormap = 'plasma';
        
        // Interaction
        this.mouseX = -1;
        this.mouseY = -1;
        
        // Setup event listeners
        this.setupEventListeners();
    }
    
    /**
     * Setup mouse interaction
     */
    setupEventListeners() {
        this.canvas.addEventListener('mousemove', (e) => {
            const rect = this.canvas.getBoundingClientRect();
            this.mouseX = e.clientX - rect.left;
            this.mouseY = e.clientY - rect.top;
            this.drawTooltip();
        });
        
        this.canvas.addEventListener('mouseleave', () => {
            this.mouseX = -1;
            this.mouseY = -1;
            this.render();
        });
    }
    
    /**
     * Set data for visualization
     */
    setData(results) {
        this.xtDataFull = results.xtData;  // Store full resolution data
        this.tracers = results.tracers;
        
        // Determine data ranges
        this.xMin = 0;
        this.xMax = results.x[results.x.length - 1] + (results.x[1] - results.x[0]);
        this.tMin = 0;
        this.tMax = this.xtDataFull[this.xtDataFull.length - 1].t;
        
        // Find pressure range from full data
        this.pMin = Infinity;
        this.pMax = -Infinity;
        for (let snapshot of this.xtDataFull) {
            for (let p of snapshot.p) {
                this.pMin = Math.min(this.pMin, p);
                this.pMax = Math.max(this.pMax, p);
            }
        }
        
        // Downsample data to canvas resolution for rendering
        this.xtData = this.downsampleData(this.xtDataFull);
    }
    
    /**
     * Downsample X-T data to match canvas resolution
     * Uses averaging to preserve data accuracy
     */
    downsampleData(fullData) {
        const fullNx = fullData[0].p.length;
        const fullNt = fullData.length;
        
        // Target resolution matches canvas pixels
        const targetNx = Math.min(this.plotWidth, fullNx);
        const targetNt = Math.min(this.plotHeight, fullNt);
        
        // If already at or below target resolution, return original data
        if (fullNx <= targetNx && fullNt <= targetNt) {
            console.log(`No downsampling needed: ${fullNx}×${fullNt} fits in ${targetNx}×${targetNt}`);
            return fullData;
        }
        
        console.log(`Downsampling from ${fullNx}×${fullNt} to ${targetNx}×${targetNt}`);
        
        const downsampled = [];
        
        // Calculate sampling factors
        const spatialFactor = fullNx / targetNx;
        const temporalFactor = fullNt / targetNt;
        
        // Downsample in time
        for (let jTarget = 0; jTarget < targetNt; jTarget++) {
            // Find source time indices to average
            const jStart = Math.floor(jTarget * temporalFactor);
            const jEnd = Math.min(Math.floor((jTarget + 1) * temporalFactor), fullNt);
            const numTimeSteps = jEnd - jStart;
            
            // Create downsampled snapshot
            const snapshot = {
                t: 0,
                p: new Float64Array(targetNx)
            };
            
            // Average time values
            for (let j = jStart; j < jEnd; j++) {
                snapshot.t += fullData[j].t;
            }
            snapshot.t /= numTimeSteps;
            
            // Downsample in space
            for (let iTarget = 0; iTarget < targetNx; iTarget++) {
                // Find source spatial indices to average
                const iStart = Math.floor(iTarget * spatialFactor);
                const iEnd = Math.min(Math.floor((iTarget + 1) * spatialFactor), fullNx);
                const numSpatialPoints = iEnd - iStart;
                
                // Average pressure over the spatial-temporal block
                let pSum = 0;
                let count = 0;
                for (let j = jStart; j < jEnd; j++) {
                    for (let i = iStart; i < iEnd; i++) {
                        pSum += fullData[j].p[i];
                        count++;
                    }
                }
                snapshot.p[iTarget] = pSum / count;
            }
            
            downsampled.push(snapshot);
        }
        
        console.log(`Downsampling complete: ${downsampled.length} snapshots with ${downsampled[0].p.length} points each`);
        return downsampled;
    }
    
    /**
     * Render the complete visualization
     */
    render() {
        // Clear canvas
        this.ctx.fillStyle = '#ffffff';
        this.ctx.fillRect(0, 0, this.width, this.height);
        
        if (!this.xtData) return;
        
        // Draw heatmap
        this.drawHeatmap();
        
        // Draw tracer lines
        this.drawTracers();
        
        // Draw axes and labels
        this.drawAxes();
        
        // Draw colorbar
        this.drawColorbar();
    }
    
    /**
     * Draw pressure heatmap
     */
    drawHeatmap() {
        const nx = this.xtData[0].p.length;
        const nt = this.xtData.length;
        
        const dx = this.plotWidth / nx;
        const dt = this.plotHeight / nt;
        
        // Create image data for faster rendering
        const imageData = this.ctx.createImageData(this.plotWidth, this.plotHeight);
        const data = imageData.data;
        
        for (let j = 0; j < nt; j++) {
            const snapshot = this.xtData[j];
            
            for (let i = 0; i < nx; i++) {
                const p = snapshot.p[i];
                
                // Normalize pressure
                let pNorm;
                if (this.useLogScale) {
                    pNorm = (Math.log10(p) - Math.log10(this.pMin)) / 
                            (Math.log10(this.pMax) - Math.log10(this.pMin));
                } else {
                    pNorm = (p - this.pMin) / (this.pMax - this.pMin);
                }
                pNorm = Math.max(0, Math.min(1, pNorm));
                
                // Get color
                const color = this.infernoColormap(pNorm);
                
                // Fill pixels (inverted y-axis: t increases upward)
                const x0 = Math.floor(i * dx);
                const x1 = Math.floor((i + 1) * dx);
                const y0 = Math.floor(this.plotHeight - (j + 1) * dt);
                const y1 = Math.floor(this.plotHeight - j * dt);
                
                for (let py = y0; py < y1; py++) {
                    for (let px = x0; px < x1; px++) {
                        const idx = (py * this.plotWidth + px) * 4;
                        data[idx] = color[0];
                        data[idx + 1] = color[1];
                        data[idx + 2] = color[2];
                        data[idx + 3] = 255;
                    }
                }
            }
        }
        
        // Draw the image
        this.ctx.putImageData(imageData, this.margin.left, this.margin.top);
    }
    
    /**
     * Draw interface tracer lines
     */
    drawTracers() {
        if (!this.tracers || this.tracers.length === 0) return;
        
        this.ctx.strokeStyle = '#ffffff';
        this.ctx.lineWidth = 2;
        this.ctx.setLineDash([]);
        
        for (let tracer of this.tracers) {
            this.ctx.beginPath();
            
            for (let i = 0; i < tracer.trajectory.length; i++) {
                const point = tracer.trajectory[i];
                const x = this.margin.left + (point.x - this.xMin) / (this.xMax - this.xMin) * this.plotWidth;
                // Inverted y-axis: t increases upward
                const y = this.margin.top + this.plotHeight - (point.t - this.tMin) / (this.tMax - this.tMin) * this.plotHeight;
                
                if (i === 0) {
                    this.ctx.moveTo(x, y);
                } else {
                    this.ctx.lineTo(x, y);
                }
            }
            
            this.ctx.stroke();
        }
    }
    
    /**
     * Draw axes and labels
     */
    drawAxes() {
        this.ctx.strokeStyle = '#000000';
        this.ctx.lineWidth = 1.5;
        this.ctx.fillStyle = '#000000';
        this.ctx.font = '14px Arial';
        this.ctx.textAlign = 'center';
        
        // Draw border
        this.ctx.strokeRect(this.margin.left, this.margin.top, this.plotWidth, this.plotHeight);
        
        // X-axis label
        this.ctx.fillText('Position (m)', this.margin.left + this.plotWidth / 2, this.height - 15);
        
        // Y-axis label
        this.ctx.save();
        this.ctx.translate(15, this.margin.top + this.plotHeight / 2);
        this.ctx.rotate(-Math.PI / 2);
        this.ctx.fillText('Time (ms)', 0, 0);
        this.ctx.restore();
        
        // Title
        this.ctx.font = 'bold 16px Arial';
        this.ctx.fillText('X-T Diagram: Pressure Field', this.width / 2, 20);
        
        // X-axis ticks
        const nTicksX = 5;
        this.ctx.font = '12px Arial';
        for (let i = 0; i <= nTicksX; i++) {
            const x = this.margin.left + i * this.plotWidth / nTicksX;
            const xVal = this.xMin + i * (this.xMax - this.xMin) / nTicksX;
            
            this.ctx.beginPath();
            this.ctx.moveTo(x, this.margin.top + this.plotHeight);
            this.ctx.lineTo(x, this.margin.top + this.plotHeight + 5);
            this.ctx.stroke();
            
            this.ctx.fillText(xVal.toFixed(2), x, this.margin.top + this.plotHeight + 20);
        }
        
        // Y-axis ticks (inverted: t increases upward)
        const nTicksY = 5;
        this.ctx.textAlign = 'right';
        for (let i = 0; i <= nTicksY; i++) {
            const y = this.margin.top + this.plotHeight - i * this.plotHeight / nTicksY;
            const tVal = this.tMin + i * (this.tMax - this.tMin) / nTicksY;
            
            this.ctx.beginPath();
            this.ctx.moveTo(this.margin.left - 5, y);
            this.ctx.lineTo(this.margin.left, y);
            this.ctx.stroke();
            
            this.ctx.fillText((tVal * 1000).toFixed(1), this.margin.left - 10, y + 4);
        }
    }
    
    /**
     * Draw colorbar
     */
    drawColorbar() {
        const barWidth = 20;
        const barHeight = this.plotHeight;
        const barX = this.width - this.margin.right + 20;
        const barY = this.margin.top;
        
        // Draw color gradient
        const gradient = this.ctx.createLinearGradient(0, barY + barHeight, 0, barY);
        for (let i = 0; i <= 10; i++) {
            const color = this.infernoColormap(i / 10);
            gradient.addColorStop(i / 10, `rgb(${color[0]},${color[1]},${color[2]})`);
        }
        
        this.ctx.fillStyle = gradient;
        this.ctx.fillRect(barX, barY, barWidth, barHeight);
        
        // Draw border
        this.ctx.strokeStyle = '#000000';
        this.ctx.lineWidth = 1;
        this.ctx.strokeRect(barX, barY, barWidth, barHeight);
        
        // Draw ticks and labels
        this.ctx.font = '11px Arial';
        this.ctx.textAlign = 'left';
        this.ctx.fillStyle = '#000000';
        
        const nTicks = 5;
        for (let i = 0; i <= nTicks; i++) {
            const y = barY + barHeight - i * barHeight / nTicks;
            const pVal = this.pMin + i * (this.pMax - this.pMin) / nTicks;
            
            // Tick mark
            this.ctx.beginPath();
            this.ctx.moveTo(barX + barWidth, y);
            this.ctx.lineTo(barX + barWidth + 5, y);
            this.ctx.stroke();
            
            // Label in kPa
            this.ctx.fillText((pVal / 1000).toFixed(0), barX + barWidth + 10, y + 4);
        }
        
        // Colorbar label
        this.ctx.save();
        this.ctx.translate(barX + barWidth + 55, barY + barHeight / 2);
        this.ctx.rotate(-Math.PI / 2);
        this.ctx.textAlign = 'center';
        this.ctx.font = '12px Arial';
        this.ctx.fillText('Pressure (kPa)', 0, 0);
        this.ctx.restore();
    }
    
    /**
     * Draw tooltip on hover
     */
    drawTooltip() {
        if (this.mouseX < 0 || this.mouseY < 0) return;
        
        // Check if mouse is in plot area
        if (this.mouseX < this.margin.left || this.mouseX > this.margin.left + this.plotWidth ||
            this.mouseY < this.margin.top || this.mouseY > this.margin.top + this.plotHeight) {
            this.render();
            return;
        }
        
        // Re-render to clear previous tooltip
        this.render();
        
        // Get data coordinates (inverted y-axis)
        const x = this.xMin + (this.mouseX - this.margin.left) / this.plotWidth * (this.xMax - this.xMin);
        const t = this.tMax - (this.mouseY - this.margin.top) / this.plotHeight * (this.tMax - this.tMin);
        
        // Find nearest data point in downsampled data
        const nx = this.xtData[0].p.length;
        const nt = this.xtData.length;
        
        const iTime = Math.max(0, Math.min(nt - 1, Math.floor((t - this.tMin) / (this.tMax - this.tMin) * nt)));
        const iPos = Math.max(0, Math.min(nx - 1, Math.floor((x - this.xMin) / (this.xMax - this.xMin) * nx)));
        
        const p = this.xtData[iTime].p[iPos];
        
        // Draw tooltip
        const tooltipText = `x: ${x.toFixed(3)} m, t: ${(t * 1000).toFixed(2)} ms, p: ${(p / 1000).toFixed(1)} kPa`;
        
        this.ctx.font = '12px Arial';
        const textWidth = this.ctx.measureText(tooltipText).width;
        
        let tooltipX = this.mouseX + 10;
        let tooltipY = this.mouseY - 10;
        
        // Keep tooltip in bounds
        if (tooltipX + textWidth + 10 > this.width) {
            tooltipX = this.mouseX - textWidth - 10;
        }
        if (tooltipY < 20) {
            tooltipY = this.mouseY + 20;
        }
        
        // Draw background
        this.ctx.fillStyle = 'rgba(255, 255, 255, 0.9)';
        this.ctx.fillRect(tooltipX - 5, tooltipY - 15, textWidth + 10, 20);
        this.ctx.strokeStyle = '#000000';
        this.ctx.lineWidth = 1;
        this.ctx.strokeRect(tooltipX - 5, tooltipY - 15, textWidth + 10, 20);
        
        // Draw text
        this.ctx.fillStyle = '#000000';
        this.ctx.textAlign = 'left';
        this.ctx.fillText(tooltipText, tooltipX, tooltipY);
        
        // Draw crosshair
        this.ctx.strokeStyle = 'rgba(255, 255, 255, 0.5)';
        this.ctx.lineWidth = 1;
        this.ctx.setLineDash([5, 5]);
        
        this.ctx.beginPath();
        this.ctx.moveTo(this.margin.left, this.mouseY);
        this.ctx.lineTo(this.margin.left + this.plotWidth, this.mouseY);
        this.ctx.stroke();
        
        this.ctx.beginPath();
        this.ctx.moveTo(this.mouseX, this.margin.top);
        this.ctx.lineTo(this.mouseX, this.margin.top + this.plotHeight);
        this.ctx.stroke();
        
        this.ctx.setLineDash([]);
    }
    
    /**
     * Inferno colormap (perceptually uniform)
     * Actual matplotlib inferno colormap data
     */
    infernoColormap(t) {
        t = Math.max(0, Math.min(1, t));
        
        // Inferno colormap lookup table (256 entries from matplotlib)
        const inferno = [
            [0,0,4],[1,0,5],[1,1,6],[1,1,8],[2,1,10],[2,2,12],[2,2,14],[3,2,16],[4,3,18],[4,3,20],[5,4,23],[6,4,25],[7,5,27],[8,5,29],[9,6,31],[10,7,34],[11,7,36],[12,8,38],[13,8,41],[14,9,43],[16,9,45],[17,10,48],[18,10,50],[20,11,52],[21,11,55],[22,11,57],[24,12,60],[25,12,62],[27,12,65],[28,12,67],[30,12,69],[31,12,72],[33,12,74],[35,12,76],[36,12,79],[38,12,81],[40,11,83],[41,11,85],[43,11,87],[45,11,89],[47,10,91],[49,10,92],[50,10,94],[52,10,95],[54,9,97],[56,9,98],[57,9,99],[59,9,100],[61,9,101],[62,9,102],[64,10,103],[66,10,104],[68,10,104],[69,10,105],[71,11,106],[73,11,106],[74,12,107],[76,12,107],[77,13,108],[79,13,108],[81,14,108],[82,14,109],[84,15,109],[85,15,109],[87,16,110],[89,16,110],[90,17,110],[92,18,110],[93,18,110],[95,19,110],[96,19,110],[98,20,110],[100,21,110],[101,21,110],[103,22,110],[104,22,110],[106,23,110],[108,24,110],[109,24,110],[111,25,110],[112,25,110],[114,26,110],[115,26,110],[117,27,110],[119,28,109],[120,28,109],[122,29,109],[123,29,109],[125,30,109],[126,30,109],[128,31,108],[129,31,108],[131,32,108],[133,32,107],[134,33,107],[136,33,107],[137,34,106],[139,34,106],[140,35,105],[142,35,105],[144,36,105],[145,36,104],[147,37,104],[148,37,103],[150,38,103],[151,38,102],[153,39,102],[154,39,101],[156,40,100],[157,40,100],[159,41,99],[160,41,99],[162,42,98],[163,42,97],[165,43,97],[166,43,96],[168,44,95],[169,44,95],[171,45,94],[172,46,93],[174,46,92],[175,47,92],[177,47,91],[178,48,90],[180,48,89],[181,49,88],[183,49,88],[184,50,87],[186,50,86],[187,51,85],[188,51,84],[190,52,83],[191,52,82],[193,53,81],[194,53,80],[196,54,80],[197,54,79],[199,55,78],[200,55,77],[201,56,76],[203,56,75],[204,57,74],[206,57,73],[207,58,72],[208,58,71],[210,59,70],[211,59,69],[213,60,68],[214,60,67],[215,61,66],[217,61,65],[218,62,64],[219,62,62],[221,63,61],[222,63,60],[223,64,59],[225,64,58],[226,65,57],[227,65,56],[228,66,55],[230,66,54],[231,67,53],[232,67,52],[234,68,51],[235,68,50],[236,69,49],[237,69,48],[239,70,47],[240,70,45],[241,71,44],[243,71,43],[244,72,42],[245,72,41],[246,73,40],[248,73,39],[249,74,38],[250,74,37],[251,75,36],[253,75,35],[254,76,34],[255,76,33],[255,77,33],[255,78,33],[255,79,33],[255,80,33],[255,81,34],[255,82,35],[255,83,35],[255,84,36],[255,85,37],[255,86,38],[255,87,39],[255,88,40],[255,89,41],[255,90,42],[255,91,43],[255,93,45],[255,94,46],[255,95,47],[255,96,48],[255,97,50],[255,98,51],[255,99,53],[255,100,54],[255,101,56],[255,102,57],[255,103,59],[255,104,61],[255,105,62],[255,106,64],[255,107,66],[255,108,68],[255,109,70],[255,110,71],[255,112,73],[255,113,75],[255,114,77],[255,115,79],[255,116,81],[255,117,83],[255,118,85],[255,119,87],[255,120,89],[255,122,91],[255,123,94],[255,124,96],[255,125,98],[255,126,100],[255,127,103],[255,128,105],[255,129,107],[255,130,109],[255,131,112],[255,132,114],[255,133,116],[255,134,119],[255,135,121],[255,136,123],[255,137,126],[255,138,128],[255,140,131],[255,141,133],[255,142,136],[255,143,138],[255,144,141],[255,145,143],[255,146,146],[255,147,148],[255,148,151],[255,149,153],[255,150,156],[255,151,159],[255,152,161],[255,153,164],[255,154,167],[255,155,169],[255,156,172],[255,157,175],[255,158,178],[255,159,180],[255,160,183],[255,161,186],[255,162,189],[255,163,192],[255,164,194],[255,165,197],[255,166,200],[255,167,203],[255,168,206],[255,169,209],[255,170,212],[255,171,215],[255,172,218],[255,173,220],[255,174,223],[255,175,226],[255,176,229],[255,177,232],[255,178,235],[255,179,238],[255,180,241],[255,181,244],[255,182,247],[255,182,250],[255,183,253],[255,184,255]
        ];
        
        // Get index in lookup table
        const idx = Math.floor(t * (inferno.length - 1));
        const color = inferno[idx];
        
        return [color[0], color[1], color[2]];
    }
    
    /**
     * Toggle log scale
     */
    toggleLogScale() {
        this.useLogScale = !this.useLogScale;
        this.render();
    }
    
    /**
     * Export canvas as PNG
     */
    exportImage(filename = 'xt-diagram.png') {
        const link = document.createElement('a');
        link.download = filename;
        link.href = this.canvas.toDataURL('image/png');
        link.click();
    }
}
