// Seeded sample spaces and classic minimal triangulations with known homology.
// Everything here is deterministic so tests and benchmarks are reproducible.

import { closure } from './filtration.js';

/** mulberry32: tiny, fast, good enough for test data. */
export function rng(seed = 1) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function gaussian(rand) {
  let u = 0;
  while (u === 0) u = rand();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rand());
}

export function circle(n, { radius = 1, noise = 0, seed = 1, center = [0, 0] } = {}) {
  const rand = rng(seed);
  return Array.from({ length: n }, (_, i) => {
    const t = (2 * Math.PI * i) / n;
    return [
      center[0] + radius * Math.cos(t) + noise * gaussian(rand),
      center[1] + radius * Math.sin(t) + noise * gaussian(rand),
    ];
  });
}

export function twoCircles(n, { noise = 0, seed = 1 } = {}) {
  const half = Math.floor(n / 2);
  return [
    ...circle(half, { radius: 1, noise, seed, center: [-1.5, 0] }),
    ...circle(n - half, { radius: 0.6, noise, seed: seed + 1, center: [1.4, 0] }),
  ];
}

export function figureEight(n, { noise = 0, seed = 1 } = {}) {
  const rand = rng(seed);
  return Array.from({ length: n }, (_, i) => {
    const t = (2 * Math.PI * i) / n;
    // lemniscate of Gerono: two lobes touching at the origin
    return [Math.sin(t) + noise * gaussian(rand), Math.sin(t) * Math.cos(t) + noise * gaussian(rand)];
  });
}

export function sphere(n, { radius = 1, seed = 1 } = {}) {
  const rand = rng(seed);
  return Array.from({ length: n }, () => {
    const v = [gaussian(rand), gaussian(rand), gaussian(rand)];
    const r = Math.hypot(...v);
    return v.map((x) => (radius * x) / r);
  });
}

export function torus(n, { R = 2, r = 1, seed = 1 } = {}) {
  const rand = rng(seed);
  return Array.from({ length: n }, () => {
    const u = 2 * Math.PI * rand();
    const v = 2 * Math.PI * rand();
    return [(R + r * Math.cos(v)) * Math.cos(u), (R + r * Math.cos(v)) * Math.sin(u), r * Math.sin(v)];
  });
}

export function uniformBox(n, dim = 2, { seed = 1 } = {}) {
  const rand = rng(seed);
  return Array.from({ length: n }, () => Array.from({ length: dim }, rand));
}

// Minimal triangulations. Over Z/2 their Betti numbers are:
//   sphere S^2 (boundary of tetrahedron)   1 0 1
//   torus T^2 (Csaszar, 7 vertices)        1 2 1
//   projective plane RP^2 (6 vertices)     1 1 1
//   Klein bottle (9 vertex grid quotient)  1 2 1
export const triangulations = {
  sphere: () => closure([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]),
  torus: () => {
    const tris = [];
    for (let i = 0; i < 7; i++) {
      tris.push([i, (i + 1) % 7, (i + 3) % 7]);
      tris.push([i, (i + 2) % 7, (i + 3) % 7]);
    }
    return closure(tris);
  },
  projectivePlane: () =>
    closure([
      [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 1],
      [1, 2, 4], [2, 3, 5], [3, 4, 1], [4, 5, 2], [5, 1, 3],
    ]),
  kleinBottle: () => gridQuotient(true),
  gridTorus: () => gridQuotient(false),
};

// 3x3 grid of squares with opposite sides identified; flipping one pair of
// sides gives the Klein bottle instead of the torus.
function gridQuotient(flip) {
  const id = (i, j) => {
    let x = ((i % 3) + 3) % 3;
    let y = ((j % 3) + 3) % 3;
    if (flip && Math.floor(i / 3) % 2 !== 0) y = (3 - y) % 3;
    return 3 * x + y;
  };
  const tris = [];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      const a = id(i, j);
      const b = id(i + 1, j);
      const c = id(i, j + 1);
      const d = id(i + 1, j + 1);
      tris.push([a, b, d], [a, c, d]);
    }
  }
  return closure(tris);
}
