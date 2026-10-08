// Distances between persistence diagrams.
//
// Both distances are optimal partial matchings where any point may instead be
// sent to its nearest point on the diagonal. The usual trick turns this into
// a square assignment problem: the left side is A plus |B| diagonal slots,
// the right side is B plus |A| diagonal slots, every point may use any
// diagonal slot at cost persistence/2 and diagonal-to-diagonal is free.
// Points at infinity are matched separately: their counts must agree, and
// sorting births is optimal for any convex cost in one dimension.

import { toPairs } from './reduce.js';

const linf = (a, b) => Math.max(Math.abs(a[0] - b[0]), Math.abs(a[1] - b[1]));
const toDiagonal = (a) => (a[1] - a[0]) / 2;

function normalise(dgm) {
  const finite = [];
  const infinite = [];
  for (const [b, d] of toPairs(dgm)) {
    if (Number.isNaN(b) || Number.isNaN(d)) throw new RangeError('diagram contains NaN');
    if (d < b) throw new RangeError(`point (${b}, ${d}) lies below the diagonal`);
    if (d === Infinity) infinite.push(b);
    else if (d > b) finite.push([b, d]);
  }
  infinite.sort((x, y) => x - y);
  return { finite, infinite };
}

export function costMatrix(A, B) {
  const n = A.length;
  const m = B.length;
  const N = n + m;
  const C = Array.from({ length: N }, () => new Float64Array(N));
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      if (i < n && j < m) C[i][j] = linf(A[i], B[j]);
      else if (i < n) C[i][j] = toDiagonal(A[i]);
      else if (j < m) C[i][j] = toDiagonal(B[j]);
      else C[i][j] = 0;
    }
  }
  return C;
}

// Hopcroft-Karp on the threshold graph { (i, j) : C[i][j] <= eps }.
function maximumMatching(C, eps) {
  const N = C.length;
  const adj = C.map((row) => {
    const out = [];
    for (let j = 0; j < N; j++) if (row[j] <= eps) out.push(j);
    return out;
  });
  const matchL = new Int32Array(N).fill(-1);
  const matchR = new Int32Array(N).fill(-1);
  const dist = new Int32Array(N);
  const INF = 1 << 30;

  const bfs = () => {
    const queue = [];
    let found = false;
    for (let u = 0; u < N; u++) {
      if (matchL[u] === -1) {
        dist[u] = 0;
        queue.push(u);
      } else dist[u] = INF;
    }
    for (let q = 0; q < queue.length; q++) {
      const u = queue[q];
      for (const v of adj[u]) {
        const w = matchR[v];
        if (w === -1) found = true;
        else if (dist[w] === INF) {
          dist[w] = dist[u] + 1;
          queue.push(w);
        }
      }
    }
    return found;
  };
  const dfs = (u) => {
    for (const v of adj[u]) {
      const w = matchR[v];
      if (w === -1 || (dist[w] === dist[u] + 1 && dfs(w))) {
        matchL[u] = v;
        matchR[v] = u;
        return true;
      }
    }
    dist[u] = INF;
    return false;
  };

  let size = 0;
  while (bfs()) {
    for (let u = 0; u < N; u++) if (matchL[u] === -1 && dfs(u)) size++;
  }
  return { size, matchL };
}

function describeMatching(matchL, n, m) {
  const out = [];
  for (let i = 0; i < matchL.length; i++) {
    const j = matchL[i];
    if (i < n) out.push([i, j < m ? j : -1]);
    else if (j < m) out.push([-1, j]);
  }
  return out;
}

/**
 * Exact bottleneck distance. The optimum is always one of the O(N^2) entries
 * of the cost matrix, so we binary search over the sorted distinct entries
 * and test each candidate for a perfect matching.
 *
 * Returns { distance, matching } where matching lists [indexInA, indexInB]
 * for the finite points (-1 means "matched to the diagonal").
 */
export function bottleneckMatching(dgmA, dgmB) {
  const A = normalise(dgmA);
  const B = normalise(dgmB);
  if (A.infinite.length !== B.infinite.length) return { distance: Infinity, matching: null };
  let infPart = 0;
  for (let k = 0; k < A.infinite.length; k++) {
    infPart = Math.max(infPart, Math.abs(A.infinite[k] - B.infinite[k]));
  }
  const n = A.finite.length;
  const m = B.finite.length;
  if (n + m === 0) return { distance: infPart, matching: [] };

  const C = costMatrix(A.finite, B.finite);
  const candidates = [...new Set(C.flatMap((row) => [...row]))].sort((x, y) => x - y);
  let lo = 0;
  let hi = candidates.length - 1;
  let best = null;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    const result = maximumMatching(C, candidates[mid]);
    if (result.size === C.length) {
      best = { eps: candidates[mid], matchL: result.matchL };
      hi = mid - 1;
    } else lo = mid + 1;
  }
  return {
    distance: Math.max(best.eps, infPart),
    matching: describeMatching(best.matchL, n, m),
    finite: { A: A.finite, B: B.finite },
  };
}

export function bottleneck(dgmA, dgmB) {
  return bottleneckMatching(dgmA, dgmB).distance;
}

/**
 * Hungarian algorithm (Kuhn-Munkres with potentials, O(N^3)) for a square
 * cost matrix. Returns { cost, assignment } with assignment[row] = column.
 */
export function hungarian(cost) {
  const n = cost.length;
  const u = new Float64Array(n + 1);
  const v = new Float64Array(n + 1);
  const p = new Int32Array(n + 1);
  const way = new Int32Array(n + 1);
  for (let i = 1; i <= n; i++) {
    p[0] = i;
    let j0 = 0;
    const minv = new Float64Array(n + 1).fill(Infinity);
    const used = new Uint8Array(n + 1);
    do {
      used[j0] = 1;
      const i0 = p[j0];
      let delta = Infinity;
      let j1 = 0;
      for (let j = 1; j <= n; j++) {
        if (used[j]) continue;
        const cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
        if (cur < minv[j]) {
          minv[j] = cur;
          way[j] = j0;
        }
        if (minv[j] < delta) {
          delta = minv[j];
          j1 = j;
        }
      }
      for (let j = 0; j <= n; j++) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        } else minv[j] -= delta;
      }
      j0 = j1;
    } while (p[j0] !== 0);
    do {
      const j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0);
  }
  const assignment = new Int32Array(n);
  let total = 0;
  for (let j = 1; j <= n; j++) {
    assignment[p[j] - 1] = j - 1;
    total += cost[p[j] - 1][j - 1];
  }
  return { cost: total, assignment };
}

/**
 * p-Wasserstein distance with the L-infinity ground metric (the convention of
 * GUDHI's and persim's defaults). p = Infinity falls back to bottleneck.
 */
export function wasserstein(dgmA, dgmB, { p = 2 } = {}) {
  if (p === Infinity) return bottleneck(dgmA, dgmB);
  if (!(p >= 1)) throw new RangeError('p must be >= 1');
  const A = normalise(dgmA);
  const B = normalise(dgmB);
  if (A.infinite.length !== B.infinite.length) return Infinity;
  let total = 0;
  for (let k = 0; k < A.infinite.length; k++) total += Math.abs(A.infinite[k] - B.infinite[k]) ** p;
  if (A.finite.length + B.finite.length > 0) {
    const C = costMatrix(A.finite, B.finite).map((row) => row.map((c) => c ** p));
    total += hungarian(C).cost;
  }
  return total ** (1 / p);
}
