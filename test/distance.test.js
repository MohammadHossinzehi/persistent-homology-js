import { test } from 'node:test';
import assert from 'node:assert/strict';
import { bottleneck, bottleneckMatching, wasserstein, hungarian, costMatrix, rips, datasets } from '../src/index.js';

const { rng, circle, uniformBox } = datasets;
const close = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) <= eps, `${a} != ${b}`);

function randomDiagram(rand, n) {
  return Array.from({ length: n }, () => {
    const b = rand() * 2;
    return [b, b + rand() * 1.5];
  });
}

function permutations(n) {
  if (n === 0) return [[]];
  const out = [];
  for (const p of permutations(n - 1)) for (let i = 0; i <= p.length; i++) out.push([...p.slice(0, i), n - 1, ...p.slice(i)]);
  return out;
}

// Reference implementations by exhaustive search over all assignments of the
// augmented cost matrix. Only feasible for tiny diagrams, which is the point.
function bruteForce(A, B, reduce) {
  const C = costMatrix(A, B);
  let best = Infinity;
  for (const perm of permutations(C.length)) best = Math.min(best, reduce(perm.map((j, i) => C[i][j])));
  return best;
}

test('bottleneck on hand checked examples', () => {
  assert.equal(bottleneck([], []), 0);
  close(bottleneck([[0, 2]], []), 1); // send to diagonal
  close(bottleneck([[0, 2]], [[0, 2.5]]), 0.5);
  close(bottleneck([[0, 2]], [[1.6, 2.4]]), 1); // cheaper to kill both than to match: max(1, 0.4)
  close(bottleneck([[0, 4]], [[1.5, 2.5]]), 1.5); // match: max(1.5, 1.5) beats killing (2)
  close(bottleneck([[0, Infinity]], [[0.3, Infinity]]), 0.3);
  assert.equal(bottleneck([[0, Infinity]], []), Infinity);
  // zero length points are on the diagonal and cost nothing
  close(bottleneck([[1, 1], [0, 1]], [[0, 1]]), 0);
});

test('Hungarian matches brute force on random square matrices', () => {
  const rand = rng(42);
  for (let trial = 0; trial < 40; trial++) {
    const n = 1 + Math.floor(rand() * 6);
    const C = Array.from({ length: n }, () => Array.from({ length: n }, () => Math.round(rand() * 100) / 10));
    let best = Infinity;
    for (const perm of permutations(n)) best = Math.min(best, perm.reduce((s, j, i) => s + C[i][j], 0));
    const { cost, assignment } = hungarian(C);
    close(cost, best);
    assert.equal(new Set(assignment).size, n);
  }
});

test('bottleneck and Wasserstein agree with exhaustive search', () => {
  const rand = rng(7);
  for (let trial = 0; trial < 60; trial++) {
    const A = randomDiagram(rand, Math.floor(rand() * 4));
    const B = randomDiagram(rand, Math.floor(rand() * 4));
    close(bottleneck(A, B), bruteForce(A, B, (c) => Math.max(0, ...c)));
    close(wasserstein(A, B, { p: 1 }), bruteForce(A, B, (c) => c.reduce((s, x) => s + x, 0)));
    close(wasserstein(A, B, { p: 2 }), Math.sqrt(bruteForce(A, B, (c) => c.reduce((s, x) => s + x * x, 0))));
  }
});

test('distances are metrics: identity, symmetry, triangle inequality', () => {
  const rand = rng(11);
  for (let trial = 0; trial < 30; trial++) {
    const [A, B, C] = [0, 1, 2].map(() => randomDiagram(rand, 2 + Math.floor(rand() * 6)));
    for (const d of [bottleneck, (x, y) => wasserstein(x, y, { p: 1 }), (x, y) => wasserstein(x, y, { p: 2 })]) {
      close(d(A, A), 0);
      close(d(A, B), d(B, A));
      assert.ok(d(A, C) <= d(A, B) + d(B, C) + 1e-9);
    }
    assert.ok(bottleneck(A, B) <= wasserstein(A, B, { p: 2 }) + 1e-9);
    assert.ok(wasserstein(A, B, { p: 2 }) <= wasserstein(A, B, { p: 1 }) + 1e-9);
  }
});

test('the matching certifies the distance', () => {
  const A = [[0, 3], [1, 1.4], [2, 5]];
  const B = [[0.2, 3.1], [2.5, 4.4]];
  const { distance, matching, finite } = bottleneckMatching(A, B);
  let worst = 0;
  for (const [i, j] of matching) {
    let c;
    if (i >= 0 && j >= 0) c = Math.max(Math.abs(finite.A[i][0] - finite.B[j][0]), Math.abs(finite.A[i][1] - finite.B[j][1]));
    else if (i >= 0) c = (finite.A[i][1] - finite.A[i][0]) / 2;
    else c = (finite.B[j][1] - finite.B[j][0]) / 2;
    worst = Math.max(worst, c);
  }
  close(worst, distance);
  assert.equal(matching.filter(([i]) => i >= 0).length, 3);
  assert.equal(matching.filter(([, j]) => j >= 0).length, 2);
});

// Stability theorem (Cohen-Steiner, Edelsbrunner, Harer; Chazal et al. for
// Rips): moving every point by at most delta moves every pairwise distance by
// at most 2 delta, so the Rips diagrams move by at most 2 delta in bottleneck.
test('stability: perturbing a point cloud by delta moves Rips diagrams by at most 2 delta', () => {
  const rand = rng(3);
  for (const delta of [0.01, 0.05, 0.2]) {
    const X = circle(30, { noise: 0.05, seed: 4 });
    const Y = X.map(([x, y]) => {
      const a = rand() * 2 * Math.PI;
      const r = rand() * delta;
      return [x + r * Math.cos(a), y + r * Math.sin(a)];
    });
    const dx = rips(X, { maxDim: 1 }).diagrams;
    const dy = rips(Y, { maxDim: 1 }).diagrams;
    for (let k = 0; k <= 1; k++) assert.ok(bottleneck(dx[k], dy[k]) <= 2 * delta + 1e-12, `H${k} delta=${delta}`);
  }
});

test('diagram distances separate a circle from a blob', () => {
  const ring = rips(circle(40, { noise: 0.03, seed: 1 }), { maxDim: 1 }).diagrams[1];
  const ring2 = rips(circle(40, { noise: 0.03, seed: 2 }), { maxDim: 1 }).diagrams[1];
  const blob = rips(uniformBox(40, 2, { seed: 1 }).map(([x, y]) => [2 * x - 1, 2 * y - 1]), { maxDim: 1 }).diagrams[1];
  assert.ok(bottleneck(ring, ring2) < bottleneck(ring, blob) / 3);
});

test('inputs are validated', () => {
  assert.throws(() => bottleneck([[2, 1]], []), /below the diagonal/);
  assert.throws(() => wasserstein([[0, 1]], [], { p: 0.5 }), /p must be/);
  close(wasserstein([[0, 2]], [], { p: Infinity }), 1);
});
