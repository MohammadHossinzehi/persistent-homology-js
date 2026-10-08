// Homology vs cohomology, with and without clearing, on Rips complexes.
// All four variants produce identical diagrams (the test suite checks this);
// the point here is how much work each one does.
// Run: node bench/clearing.js

import { euclideanDistances, ripsFiltration, buildFiltration, persistence, datasets } from '../src/index.js';

const cases = [
  ['circle n=80, H<=1', datasets.circle(80, { noise: 0.05 }), { maxDim: 1 }],
  ['uniform box n=60, H<=1', datasets.uniformBox(60, 2), { maxDim: 1 }],
  ['sphere n=70, H<=2, r<=1.2', datasets.sphere(70, { seed: 2 }), { maxDim: 2, threshold: 1.2 }],
  ['torus n=200, H<=2, r<=1.4', datasets.torus(200, { seed: 3 }), { maxDim: 2, threshold: 1.4 }],
];

const variants = [
  ['hom', { algorithm: 'homology', clearing: false }],
  ['hom+clear', { algorithm: 'homology', clearing: true }],
  ['cohom', { algorithm: 'cohomology', clearing: false }],
  ['cohom+clear', { algorithm: 'cohomology', clearing: true }],
];

const work = [];
const time = [];
for (const [name, points, opts] of cases) {
  const filtration = buildFiltration(ripsFiltration(euclideanDistances(points), opts));
  const w = { case: name, simplices: filtration.simplices.length };
  const t = { case: name };
  for (const [label, v] of variants) {
    const t0 = performance.now();
    const r = persistence(filtration, { ...v, maxDim: opts.maxDim });
    t[label] = `${Math.round(performance.now() - t0)} ms`;
    w[label] = r.columnAdditions;
  }
  work.push(w);
  time.push(t);
}
console.log('column additions');
console.table(work);
console.log('wall time');
console.table(time);
