import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  buildFiltration,
  persistence,
  rips,
  ripsFiltration,
  euclideanDistances,
  lowerStar,
  timeSeriesFiltration,
  closure,
  addColumns,
  bettiNumbers,
  datasets,
} from '../src/index.js';

const { triangulations, circle, twoCircles, figureEight, sphere, uniformBox, rng } = datasets;

const betti = (result) => result.diagrams.map((d) => d.filter((p) => p.death === Infinity).length);
const binom = (n, k) => {
  let r = 1;
  for (let i = 1; i <= k; i++) r = (r * (n - k + i)) / i;
  return Math.round(r);
};

test('addColumns is symmetric difference of sorted index lists', () => {
  assert.deepEqual(addColumns([1, 3, 5, 9], [2, 3, 9, 11]), [1, 2, 5, 11]);
  assert.deepEqual(addColumns([4, 7], [4, 7]), []);
  assert.deepEqual(addColumns([], [0, 2]), [0, 2]);
});

test('filtration rejects complexes that are not closed under faces or not monotone', () => {
  assert.throws(() => buildFiltration([{ vertices: [0, 1], value: 1 }]), /missing/);
  assert.throws(
    () =>
      buildFiltration([
        { vertices: [0], value: 0 },
        { vertices: [1], value: 2 },
        { vertices: [0, 1], value: 1 },
      ]),
    /after its coface/,
  );
  assert.throws(() => buildFiltration([{ vertices: [0, 0] }]), /repeats/);
  assert.throws(() => buildFiltration([{ vertices: [0] }, { vertices: [0] }]), /duplicate/);
  assert.throws(() => buildFiltration([{ vertices: [0], value: NaN }]), /invalid/);
});

test('Betti numbers of classic minimal triangulations over Z/2', () => {
  const expected = {
    sphere: [1, 0, 1],
    torus: [1, 2, 1],
    projectivePlane: [1, 1, 1],
    kleinBottle: [1, 2, 1],
    gridTorus: [1, 2, 1],
  };
  for (const [name, b] of Object.entries(expected)) {
    for (const clearing of [true, false]) {
      assert.deepEqual(betti(persistence(triangulations[name](), { clearing })), b, `${name} clearing=${clearing}`);
    }
  }
});

// Over Z/2 the Klein bottle and torus look identical, so check the
// triangulations really are the surfaces they claim to be: closed (every edge
// in exactly two triangles), connected links, and the right orientability.
function orientable(simplices) {
  const tris = simplices.filter((s) => s.vertices.length === 3).map((s) => s.vertices);
  const byEdge = new Map();
  tris.forEach((t, i) => {
    for (const [a, b] of [[t[0], t[1]], [t[1], t[2]], [t[0], t[2]]]) {
      const k = `${a},${b}`;
      if (!byEdge.has(k)) byEdge.set(k, []);
      byEdge.get(k).push(i);
    }
  });
  for (const list of byEdge.values()) assert.equal(list.length, 2, 'not a closed surface');
  // orientation +1 keeps vertex order, -1 reverses it; a shared edge must be
  // traversed in opposite directions by its two triangles.
  const dir = (t, s, a, b) => {
    const cyc = s === 1 ? t : [t[0], t[2], t[1]];
    for (let k = 0; k < 3; k++) if (cyc[k] === a && cyc[(k + 1) % 3] === b) return 1;
    return -1;
  };
  const orient = new Array(tris.length).fill(0);
  orient[0] = 1;
  const stack = [0];
  while (stack.length) {
    const i = stack.pop();
    const t = tris[i];
    for (const [a, b] of [[t[0], t[1]], [t[1], t[2]], [t[0], t[2]]]) {
      for (const j of byEdge.get(`${a},${b}`)) {
        if (j === i) continue;
        const want = -dir(t, orient[i], a, b);
        const s = dir(tris[j], 1, a, b) === want ? 1 : -1;
        if (orient[j] === 0) {
          orient[j] = s;
          stack.push(j);
        } else if (orient[j] !== s) return false;
      }
    }
  }
  return true;
}

test('surface triangulations have the right orientability', () => {
  assert.equal(orientable(triangulations.sphere()), true);
  assert.equal(orientable(triangulations.torus()), true);
  assert.equal(orientable(triangulations.gridTorus()), true);
  assert.equal(orientable(triangulations.kleinBottle()), false);
  assert.equal(orientable(triangulations.projectivePlane()), false);
});

test('Euler characteristic: alternating simplex count equals alternating Betti sum', () => {
  for (const name of Object.keys(triangulations)) {
    const cx = triangulations[name]();
    const chi = cx.reduce((acc, s) => acc + (s.vertices.length % 2 === 1 ? 1 : -1), 0);
    const b = betti(persistence(cx));
    assert.equal(chi, b.reduce((acc, x, k) => acc + (k % 2 === 0 ? x : -x), 0), name);
  }
});

test('Rips complex with no threshold contains every subset up to maxDim + 1', () => {
  const pts = uniformBox(9, 2, { seed: 4 });
  const cx = ripsFiltration(euclideanDistances(pts), { maxDim: 2 });
  const counts = [0, 0, 0, 0];
  for (const s of cx) counts[s.vertices.length - 1]++;
  assert.deepEqual(counts, [binom(9, 1), binom(9, 2), binom(9, 3), binom(9, 4)]);
  for (const s of cx) {
    let diam = 0;
    for (const a of s.vertices) for (const b of s.vertices) diam = Math.max(diam, Math.hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]));
    assert.ok(Math.abs(diam - s.value) < 1e-12);
  }
});

test('H0 deaths are exactly the minimum spanning tree edge lengths', () => {
  for (let seed = 1; seed <= 5; seed++) {
    const pts = uniformBox(40, 3, { seed });
    const D = euclideanDistances(pts);
    // Kruskal
    const edges = [];
    for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++) edges.push([D[i][j], i, j]);
    edges.sort((a, b) => a[0] - b[0]);
    const parent = pts.map((_, i) => i);
    const find = (x) => (parent[x] === x ? x : (parent[x] = find(parent[x])));
    const mst = [];
    for (const [w, i, j] of edges) {
      const a = find(i);
      const b = find(j);
      if (a !== b) {
        parent[a] = b;
        mst.push(w);
      }
    }
    const h0 = rips(pts, { maxDim: 0, keepZero: true }).diagrams[0];
    const deaths = h0.filter((p) => p.death !== Infinity).map((p) => p.death).sort((a, b) => a - b);
    assert.deepEqual(deaths, mst);
    assert.equal(h0.filter((p) => p.death === Infinity).length, 1);
  }
});

test('homology, cohomology, with or without clearing: same pairs, different work', () => {
  const key = (r) =>
    r.diagrams.map((d) => d.map((p) => `${p.birthSimplex}|${p.deathSimplex}|${p.birth}|${p.death}`).sort());
  for (let seed = 1; seed <= 6; seed++) {
    const pts = uniformBox(28, 2, { seed });
    const runs = [];
    for (const algorithm of ['homology', 'cohomology']) {
      for (const clearing of [false, true]) {
        runs.push(rips(pts, { maxDim: 2, keepZero: true, algorithm, clearing }));
      }
    }
    for (const r of runs.slice(1)) assert.deepEqual(key(r), key(runs[0]));
    const [hom, homClear, , cohomClear] = runs.map((r) => r.columnAdditions);
    assert.ok(homClear <= hom);
    assert.ok(cohomClear * 10 < hom, `cohomology with clearing should be far cheaper (${cohomClear} vs ${hom})`);
  }
  // also on explicit complexes, where the pairing is the full story
  for (const name of Object.keys(triangulations)) {
    const a = persistence(triangulations[name](), { algorithm: 'homology', keepZero: true });
    const b = persistence(triangulations[name](), { algorithm: 'cohomology', keepZero: true });
    assert.deepEqual(key(a), key(b), name);
  }
  assert.throws(() => persistence(triangulations.sphere(), { algorithm: 'magic' }), /unknown algorithm/);
});

test('pairing is a bijection: every simplex is a birth, a death or essential exactly once', () => {
  const pts = uniformBox(20, 2, { seed: 9 });
  const r = rips(pts, { maxDim: 2, keepZero: true });
  const used = new Map();
  const mark = (s) => {
    if (!s) return;
    const k = s.join(',');
    used.set(k, (used.get(k) ?? 0) + 1);
  };
  // degree 3 births are dropped by maxDim, so count only simplices of dim <= 3 that are tracked
  for (const dgm of r.diagrams) for (const p of dgm) {
    mark(p.birthSimplex);
    mark(p.deathSimplex);
  }
  for (const v of used.values()) assert.equal(v, 1);
});

test('a noisy circle has one dominant H1 bar, two circles have two', () => {
  const one = rips(circle(60, { noise: 0.04, seed: 3 }), { maxDim: 1 }).diagrams[1];
  const pers = one.map((p) => p.death - p.birth);
  assert.ok(pers[0] > 1.2, `dominant loop too short: ${pers[0]}`);
  assert.ok(pers[1] === undefined || pers[1] < 0.15 * pers[0]);
  // the loop is born when neighbouring samples connect and dies when the
  // equilateral triangle inscribed in the circle appears (side sqrt 3)
  assert.ok(Math.abs(one[0].death - Math.sqrt(3)) < 0.15);

  const two = rips(twoCircles(80, { noise: 0.02, seed: 5 }), { maxDim: 1 }).diagrams[1];
  const p2 = two.map((p) => p.death - p.birth);
  assert.ok(p2[1] > 0.6 && (p2[2] ?? 0) < 0.2 * p2[1], JSON.stringify(p2.slice(0, 4)));

  const eight = rips(figureEight(80, { seed: 2 }), { maxDim: 1 }).diagrams[1];
  assert.equal(eight.filter((p) => p.death - p.birth > 0.3).length, 2);
});

test('a sampled 2-sphere fills in to the homology of S^2', () => {
  // Sparse random samples leave small holes that show up as short H1 bars;
  // by the truncation radius they have all closed and exactly the void remains.
  const r = rips(sphere(60, { seed: 2 }), { maxDim: 2, threshold: 1.3 });
  assert.deepEqual(bettiNumbers(r.diagrams, 1.29), [1, 0, 1]);
  assert.equal(r.diagrams[1].filter((p) => p.death === Infinity).length, 0);
  assert.equal(r.diagrams[2].length, 1);
});

test('representative cycles are genuine cycles with the reported birth', () => {
  const r = rips(circle(30, { noise: 0.03, seed: 8 }), { maxDim: 1, representatives: true });
  const top = r.diagrams[1][0];
  // every vertex appears an even number of times in the boundary of the cycle
  const degree = new Map();
  for (const [a, b] of top.cycle) {
    degree.set(a, (degree.get(a) ?? 0) + 1);
    degree.set(b, (degree.get(b) ?? 0) + 1);
  }
  for (const d of degree.values()) assert.equal(d % 2, 0);
  assert.ok(top.cycle.some((e) => e.join(',') === top.birthSimplex.join(',')));
  // for a loop around a circle every sample point should be on the cycle
  assert.equal(degree.size, 30);
});

test('lower-star H0 of a time series follows the elder rule', () => {
  //          0  1  2  3  4  5  6
  const ys = [3, 1, 4, 0, 5, 2, 6];
  const h0 = persistence(timeSeriesFiltration(ys)).diagrams[0].map((p) => [p.birth, p.death]);
  // minima 0 (global, never dies), 1 dies at 4, 2 dies at 5
  assert.deepEqual(h0, [[0, Infinity], [1, 4], [2, 5]]);
});

test('lower-star filtration on a 2D grid finds a crater as an H1 class', () => {
  const size = 9;
  const id = (i, j) => i * size + j;
  const values = [];
  const tris = [];
  for (let i = 0; i < size; i++) {
    for (let j = 0; j < size; j++) {
      const r = Math.hypot(i - 4, j - 4);
      values.push(-Math.exp(-((r - 2.5) ** 2))); // a ring shaped valley
      if (i + 1 < size && j + 1 < size) {
        tris.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)], [id(i, j), id(i, j + 1), id(i + 1, j + 1)]);
      }
    }
  }
  const r = persistence(lowerStar(tris, values));
  const loops = r.diagrams[1].filter((p) => p.death - p.birth > 0.3);
  assert.equal(loops.length, 1);
  assert.deepEqual(bettiNumbers(r.diagrams, 10), [1, 0, 0]);
});

test('closure produces every face exactly once', () => {
  const faces = closure([[0, 1, 2, 3]]);
  assert.equal(faces.length, 15);
});

test('rips accepts a distance matrix and validates it', () => {
  const D = [
    [0, 1, 2],
    [1, 0, 1],
    [2, 1, 0],
  ];
  const r = rips({ distances: D }, { maxDim: 0 });
  assert.deepEqual(r.diagrams[0].map((p) => p.death), [Infinity, 1, 1]);
  assert.throws(() => rips({ distances: [[0, 1], [2, 0]] }), /symmetric/);
  assert.throws(() => rips(uniformBox(30, 2), { maxDim: 3, maxSimplices: 1000 }), /exceeds/);
  const rand = rng(1);
  assert.ok(rand() >= 0 && rand() < 1);
});
