// Persistent homology over Z/2 by boundary matrix reduction.
//
// Columns are stored sparsely as ascending arrays of row indices, so the
// pivot ("low") of a column is its last entry and adding two columns is a
// linear merge that drops common entries (1 + 1 = 0 in Z/2).

import { buildFiltration } from './filtration.js';

const EMPTY = Object.freeze([]);

export function addColumns(a, b) {
  const out = [];
  let i = 0;
  let j = 0;
  while (i < a.length && j < b.length) {
    if (a[i] < b[j]) out.push(a[i++]);
    else if (a[i] > b[j]) out.push(b[j++]);
    else {
      i++;
      j++;
    }
  }
  while (i < a.length) out.push(a[i++]);
  while (j < b.length) out.push(b[j++]);
  return out;
}

/**
 * Reduce the boundary matrix of a filtration.
 *
 * With clearing (Chen and Kerber 2011, a.k.a. the twist algorithm) columns
 * are processed from the top dimension down. Whenever column j reduces to a
 * pivot i, column i is known to reduce to zero, so we skip it entirely. On
 * Rips complexes this removes the bulk of the work because most edges are
 * negative (they kill a component) or are paired with a triangle.
 *
 * Returns { R, low, pairs, essential, columnAdditions } where pairs are
 * index pairs [birth, death] in the filtration order.
 */
export function reduceBoundary(filtration, { clearing = true } = {}) {
  const { boundary, dims, maxDim } = filtration;
  const n = boundary.length;
  const R = new Array(n);
  const owner = new Int32Array(n).fill(-1); // owner[row] = column whose pivot is row
  const cleared = new Uint8Array(n);
  let columnAdditions = 0;

  const reduceColumn = (j) => {
    let col = boundary[j];
    while (col.length) {
      const k = owner[col[col.length - 1]];
      if (k === -1) break;
      col = addColumns(col, R[k]);
      columnAdditions++;
    }
    R[j] = col;
    if (col.length) {
      const pivot = col[col.length - 1];
      owner[pivot] = j;
      if (clearing) cleared[pivot] = 1;
    }
  };

  if (clearing) {
    const byDim = Array.from({ length: maxDim + 1 }, () => []);
    for (let j = 0; j < n; j++) byDim[dims[j]].push(j);
    for (let d = maxDim; d >= 0; d--) {
      for (const j of byDim[d]) {
        if (cleared[j]) R[j] = [];
        else reduceColumn(j);
      }
    }
  } else {
    for (let j = 0; j < n; j++) reduceColumn(j);
  }

  const pairs = [];
  const essential = [];
  for (let j = 0; j < n; j++) {
    if (R[j].length) pairs.push([R[j][R[j].length - 1], j]);
    else if (owner[j] === -1) essential.push(j);
  }
  pairs.sort((a, b) => a[0] - b[0]);
  return { R, owner, pairs, essential, columnAdditions };
}

/**
 * Persistent cohomology: reduce the anti-transpose of the boundary matrix
 * (the coboundary matrix with rows and columns in reverse filtration order).
 * de Silva, Morozov and Vejdemo-Johansson (2011) showed it yields exactly the
 * same persistence pairs. Combined with clearing, now applied from low
 * dimension upwards, it is the main trick behind Ripser: in a Rips complex
 * almost every simplex of dimension d + 1 is negative, so clearing removes
 * almost all of the expensive high dimensional columns.
 */
export function reduceCoboundary(filtration, { clearing = true, upTo = Infinity } = {}) {
  const { boundary, dims, maxDim } = filtration;
  const top = Math.min(maxDim, upTo);
  const n = boundary.length;
  // anti-transposed column for simplex i: rows n-1-j for every coface j,
  // stored ascending (so cofaces in descending filtration order).
  const cobound = Array.from({ length: n }, () => []);
  for (let j = n - 1; j >= 0; j--) {
    for (const f of boundary[j]) cobound[f].push(n - 1 - j);
  }
  const R = Array.from({ length: n }, () => EMPTY); // indexed by anti-transposed column c = n-1-i
  const owner = new Int32Array(n).fill(-1);
  const cleared = new Uint8Array(n);
  let columnAdditions = 0;

  const byDim = Array.from({ length: maxDim + 1 }, () => []);
  for (let i = n - 1; i >= 0; i--) byDim[dims[i]].push(i); // increasing column index
  for (let d = 0; d <= top; d++) {
    for (const i of byDim[d]) {
      const c = n - 1 - i;
      if (clearing && cleared[c]) continue;
      let col = cobound[i];
      while (col.length) {
        const k = owner[col[col.length - 1]];
        if (k === -1) break;
        col = addColumns(col, R[k]);
        columnAdditions++;
      }
      R[c] = col;
      if (col.length) {
        const pivot = col[col.length - 1];
        owner[pivot] = c;
        if (clearing) cleared[pivot] = 1;
      }
    }
  }

  const pairs = [];
  const essential = [];
  for (let c = 0; c < n; c++) {
    if (R[c].length) pairs.push([n - 1 - c, n - 1 - R[c][R[c].length - 1]]);
    else if (owner[c] === -1) essential.push(n - 1 - c);
  }
  pairs.sort((a, b) => a[0] - b[0]);
  essential.sort((a, b) => a - b);
  return { pairs, essential, columnAdditions };
}

/**
 * Compute persistence diagrams.
 *
 * `input` is either a filtration from buildFiltration or a raw array of
 * { vertices, value } simplices. Options:
 *   algorithm       'cohomology' (default) or 'homology'; same output
 *   clearing        use the clearing / twist optimisation (default true)
 *   keepZero        keep pairs with birth === death (default false)
 *   maxDim          only report homology up to this degree (default: all)
 *   representatives attach a cycle (list of simplices) to every finite
 *                   pair of degree >= 1: the reduced column of the killing
 *                   simplex, which is a cycle born exactly at the birth
 *                   simplex and killed at the death simplex. Forces the
 *                   homology algorithm, since cohomology yields cocycles.
 */
export function persistence(input, options = {}) {
  const filtration = Array.isArray(input) ? buildFiltration(input) : input;
  const { clearing = true, keepZero = false, representatives = false } = options;
  const algorithm = representatives ? 'homology' : options.algorithm ?? 'cohomology';
  if (algorithm !== 'homology' && algorithm !== 'cohomology') {
    throw new RangeError(`unknown algorithm "${algorithm}"`);
  }
  const maxDim = options.maxDim ?? filtration.maxDim;
  const { simplices, dims } = filtration;
  const reduction =
    algorithm === 'homology'
      ? reduceBoundary(filtration, { clearing })
      : reduceCoboundary(filtration, { clearing, upTo: maxDim });

  const diagrams = Array.from({ length: Math.max(0, maxDim + 1) }, () => []);
  for (const [i, j] of reduction.pairs) {
    const dim = dims[i];
    if (dim > maxDim) continue;
    const birth = simplices[i].value;
    const death = simplices[j].value;
    if (!keepZero && birth === death) continue;
    const point = {
      dim,
      birth,
      death,
      birthSimplex: simplices[i].vertices,
      deathSimplex: simplices[j].vertices,
    };
    if (representatives && dim >= 1) point.cycle = reduction.R[j].map((r) => simplices[r].vertices);
    diagrams[dim].push(point);
  }
  for (const i of reduction.essential) {
    const dim = dims[i];
    if (dim > maxDim) continue;
    diagrams[dim].push({
      dim,
      birth: simplices[i].value,
      death: Infinity,
      birthSimplex: simplices[i].vertices,
      deathSimplex: null,
    });
  }
  for (const dgm of diagrams) {
    dgm.sort((a, b) => b.death - b.birth - (a.death - a.birth) || a.birth - b.birth);
  }
  return {
    diagrams,
    numSimplices: simplices.length,
    algorithm,
    columnAdditions: reduction.columnAdditions,
    filtration,
  };
}

/** Strip a diagram to plain [birth, death] pairs, the input format of the distance functions. */
export function toPairs(diagram) {
  return diagram.map((p) => (Array.isArray(p) ? [p[0], p[1]] : [p.birth, p.death]));
}
