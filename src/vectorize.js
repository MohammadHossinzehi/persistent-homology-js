// Summaries of persistence diagrams: Betti curves, persistence landscapes,
// persistence images and persistent entropy. These turn a diagram into a
// fixed length vector you can feed to any ordinary ML model.

import { toPairs } from './reduce.js';

/** Number of bars alive at time t (birth <= t < death). */
export function bettiNumber(diagram, t) {
  let count = 0;
  for (const [b, d] of toPairs(diagram)) if (b <= t && t < d) count++;
  return count;
}

/** Betti numbers for every degree of a persistence() result at time t. */
export function bettiNumbers(diagrams, t) {
  return diagrams.map((dgm) => bettiNumber(dgm, t));
}

/** Betti curve sampled on a grid. */
export function bettiCurve(diagram, grid) {
  const pairs = toPairs(diagram);
  return grid.map((t) => pairs.reduce((acc, [b, d]) => acc + (b <= t && t < d ? 1 : 0), 0));
}

export function linspace(start, stop, count) {
  if (count < 2) return [start];
  const step = (stop - start) / (count - 1);
  return Array.from({ length: count }, (_, i) => start + i * step);
}

/**
 * Persistence landscape (Bubenik 2015). Each bar (b, d) contributes the tent
 * max(0, min(t - b, d - t)); lambda_k(t) is the k-th largest tent value at t.
 * Infinite bars are clipped to `infinity` (default: last grid point), which
 * is the usual convention since landscapes must be integrable.
 *
 * Returns an array of `levels` arrays, each sampled on `grid`.
 */
export function landscape(diagram, grid, { levels = 3, infinity } = {}) {
  const cap = infinity ?? grid[grid.length - 1];
  const bars = toPairs(diagram).map(([b, d]) => [b, d === Infinity ? Math.max(cap, b) : d]);
  const out = Array.from({ length: levels }, () => new Float64Array(grid.length));
  const buf = new Float64Array(bars.length);
  for (let g = 0; g < grid.length; g++) {
    const t = grid[g];
    for (let i = 0; i < bars.length; i++) {
      const [b, d] = bars[i];
      buf[i] = Math.max(0, Math.min(t - b, d - t));
    }
    const sorted = Array.from(buf).sort((x, y) => y - x);
    for (let k = 0; k < levels; k++) out[k][g] = sorted[k] ?? 0;
  }
  return out.map((a) => Array.from(a));
}

/**
 * Persistence image (Adams et al. 2017). Points are mapped to
 * (birth, persistence), weighted by a ramp that is 0 at persistence 0 and 1 at
 * the largest persistence, then smoothed by an isotropic Gaussian and
 * integrated over pixels by evaluating the Gaussian CDF at pixel edges, so the
 * image is exact rather than point sampled.
 */
export function persistenceImage(diagram, { resolution = 20, sigma = 0.1, birthRange, persistenceRange } = {}) {
  const pts = toPairs(diagram)
    .filter(([, d]) => d !== Infinity)
    .map(([b, d]) => [b, d - b]);
  const maxPers = pts.reduce((m, [, p]) => Math.max(m, p), 0) || 1;
  const [bx0, bx1] = birthRange ?? [Math.min(0, ...pts.map((q) => q[0])), Math.max(1e-9, ...pts.map((q) => q[0]))];
  const [py0, py1] = persistenceRange ?? [0, maxPers];
  const erf = (x) => {
    // Abramowitz and Stegun 7.1.26, |error| < 1.5e-7
    const s = Math.sign(x);
    const ax = Math.abs(x);
    const t = 1 / (1 + 0.3275911 * ax);
    const y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-ax * ax);
    return s * y;
  };
  const cdf = (x, mu) => 0.5 * (1 + erf((x - mu) / (sigma * Math.SQRT2)));
  const xs = linspace(bx0, bx1 + (bx1 === bx0 ? 1 : 0), resolution + 1);
  const ys = linspace(py0, py1, resolution + 1);
  const img = Array.from({ length: resolution }, () => new Array(resolution).fill(0));
  for (const [b, p] of pts) {
    const w = p / maxPers;
    for (let r = 0; r < resolution; r++) {
      const wy = cdf(ys[r + 1], p) - cdf(ys[r], p);
      if (wy === 0) continue;
      for (let c = 0; c < resolution; c++) {
        img[r][c] += w * wy * (cdf(xs[c + 1], b) - cdf(xs[c], b));
      }
    }
  }
  return img;
}

/** Shannon entropy of the normalised bar lengths (finite bars only). */
export function persistentEntropy(diagram) {
  const lengths = toPairs(diagram)
    .filter(([b, d]) => d !== Infinity && d > b)
    .map(([b, d]) => d - b);
  const total = lengths.reduce((a, x) => a + x, 0);
  if (total === 0) return 0;
  return -lengths.reduce((acc, l) => acc + (l / total) * Math.log(l / total), 0);
}

/** Total persistence: sum of (death - birth)^p over finite bars. */
export function totalPersistence(diagram, p = 1) {
  return toPairs(diagram)
    .filter(([, d]) => d !== Infinity)
    .reduce((acc, [b, d]) => acc + (d - b) ** p, 0);
}
