/**
 * Timing read off an x-t simulation: when a wave arrives at a location and
 * how long the state behind it stays steady (the usable test time).
 *
 * Pure functions: no DOM access. Times are resolved to the snapshot interval.
 */

/**
 * First steady period in a pressure history.
 * @param {ArrayLike<number>} times
 * @param {ArrayLike<number>} values
 * @param {Object} [options]
 * @param {number} [options.arrival=0.02] - relative change from the initial
 *   value that counts as a wave arriving
 * @param {number} [options.settled=0.01] - relative change between samples
 *   below which the rise is over
 * @param {number} [options.tolerance=0.05] - relative departure from the
 *   plateau that ends the steady period
 * @returns {null|{arrival:number, start:number, plateau:number, duration:number, open:boolean}}
 *   null if nothing arrives; open is true when the plateau lasts to the end
 *   of the record, so duration is only a lower bound
 */
function steadyWindow(times, values, options) {
    const o = options || {};
    const arrivalLevel = o.arrival === undefined ? 0.02 : o.arrival;
    const settled = o.settled === undefined ? 0.01 : o.settled;
    const tolerance = o.tolerance === undefined ? 0.05 : o.tolerance;
    const n = values.length;
    const initial = values[0];

    let first = -1;
    for (let j = 1; j < n; j++) {
        if (Math.abs(values[j] - initial) > arrivalLevel * Math.abs(initial)) {
            first = j;
            break;
        }
    }
    if (first < 0) return null;

    // Step past the rise, which is spread over a few samples
    let start = first;
    while (start < n - 1 && Math.abs(values[start + 1] - values[start]) > settled * Math.abs(values[start + 1])) start++;
    const plateau = values[start];

    let end = n - 1, open = true;
    for (let j = start + 1; j < n; j++) {
        if (Math.abs(values[j] - plateau) > tolerance * Math.abs(plateau)) {
            end = j;
            open = false;
            break;
        }
    }
    return {
        arrival: times[first],
        start: times[start],
        plateau: plateau,
        duration: times[end] - times[start],
        open: open
    };
}

/** Time at which a tracer (material interface) first reaches position x, or null. */
function interfaceArrival(trajectory, x) {
    const side = trajectory[0].x < x;
    for (let j = 1; j < trajectory.length; j++) {
        if ((trajectory[j].x < x) !== side) return trajectory[j].t;
    }
    return null;
}

/**
 * Test-time summary of a simulation.
 * @param {Object} results - from EulerSolver.getResults()
 * @param {number[]} probes - positions (m)
 * @returns {{endWall: Object|null, probes: Array}} endWall is the steady
 *   period of the reflected state at the driven end wall; each probe entry is
 *   {x, window, interfaceArrival}
 */
function analyzeTestTime(results, probes) {
    const times = results.xtData.map(function (snapshot) { return snapshot.t; });
    const dx = results.x[1] - results.x[0];
    const history = function (cell) {
        return results.xtData.map(function (snapshot) { return snapshot.p[cell]; });
    };
    const last = results.x.length - 1;
    return {
        endWall: steadyWindow(times, history(last)),
        probes: (probes || []).map(function (x) {
            const cell = Math.max(0, Math.min(last, Math.floor(x / dx)));
            let arrival = null;
            for (const tracer of results.tracers || []) {
                const t = interfaceArrival(tracer.trajectory, x);
                if (t !== null && (arrival === null || t < arrival)) arrival = t;
            }
            return { x: x, window: steadyWindow(times, history(cell)), interfaceArrival: arrival };
        })
    };
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { steadyWindow, interfaceArrival, analyzeTestTime };
}
