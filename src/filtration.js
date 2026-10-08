// Filtered simplicial complexes: validation, ordering and boundary matrices,
// plus the two constructions people actually use in practice: Vietoris-Rips
// on a point cloud / distance matrix, and lower-star filtrations of a function
// defined on the vertices of a complex.

const keyOf = (vertices) => vertices.join(',');

/**
 * Total order used for the boundary matrix: by filtration value, then by
 * dimension (so a face always precedes its cofaces when they appear at the
 * same time), then lexicographically by vertex list so results are
 * deterministic across runs and platforms.
 */
export function compareSimplices(a, b) {
  if (a.value !== b.value) return a.value - b.value;
  const la = a.vertices.length;
  const lb = b.vertices.length;
  if (la !== lb) return la - lb;
  for (let i = 0; i < la; i++) {
    if (a.vertices[i] !== b.vertices[i]) return a.vertices[i] - b.vertices[i];
  }
  return 0;
}

/**
 * Build a filtration from a list of { vertices, value } records.
 *
 * Every face of every simplex must also be present, with a filtration value
 * no larger than the simplex itself; otherwise the input is not a filtered
 * simplicial complex and persistence is meaningless, so we throw.
 *
 * Returns { simplices, boundary, dims, maxDim } where boundary[j] is the
 * sorted list of indices of the codimension-1 faces of simplex j.
 */
export function buildFiltration(input) {
  if (!Array.isArray(input)) throw new TypeError('buildFiltration expects an array of simplices');
  const simplices = input.map((s, idx) => {
    if (!s || !Array.isArray(s.vertices) || s.vertices.length === 0) {
      throw new TypeError(`simplex #${idx} has no vertices`);
    }
    const value = s.value === undefined ? 0 : s.value;
    if (typeof value !== 'number' || Number.isNaN(value) || value === Infinity) {
      throw new RangeError(`simplex #${idx} has invalid filtration value ${value}`);
    }
    const vertices = [...s.vertices].sort((x, y) => x - y);
    for (let i = 0; i < vertices.length; i++) {
      if (!Number.isInteger(vertices[i]) || vertices[i] < 0) {
        throw new TypeError(`simplex #${idx} has a non integer vertex`);
      }
      if (i > 0 && vertices[i] === vertices[i - 1]) {
        throw new TypeError(`simplex #${idx} repeats vertex ${vertices[i]}`);
      }
    }
    return { vertices, value };
  });

  simplices.sort(compareSimplices);

  const index = new Map();
  for (let i = 0; i < simplices.length; i++) {
    const k = keyOf(simplices[i].vertices);
    if (index.has(k)) throw new Error(`duplicate simplex [${k}]`);
    index.set(k, i);
  }

  const boundary = new Array(simplices.length);
  const dims = new Int32Array(simplices.length);
  let maxDim = -1;
  for (let j = 0; j < simplices.length; j++) {
    const { vertices, value } = simplices[j];
    const d = vertices.length - 1;
    dims[j] = d;
    if (d > maxDim) maxDim = d;
    if (d === 0) {
      boundary[j] = [];
      continue;
    }
    const faces = new Array(d + 1);
    for (let drop = 0; drop <= d; drop++) {
      const face = vertices.filter((_, t) => t !== drop);
      const fi = index.get(keyOf(face));
      if (fi === undefined) throw new Error(`face [${face}] of simplex [${vertices}] is missing`);
      if (simplices[fi].value > value) {
        throw new Error(`face [${face}] enters at ${simplices[fi].value}, after its coface [${vertices}] at ${value}`);
      }
      faces[drop] = fi;
    }
    faces.sort((x, y) => x - y);
    boundary[j] = faces;
  }

  return { simplices, boundary, dims, maxDim };
}

/** All faces (including the simplices themselves) of a list of maximal simplices. */
export function closure(maximal, value = 0) {
  const seen = new Map();
  const visit = (verts) => {
    const k = keyOf(verts);
    if (seen.has(k)) return;
    seen.set(k, verts);
    if (verts.length > 1) {
      for (let drop = 0; drop < verts.length; drop++) visit(verts.filter((_, t) => t !== drop));
    }
  };
  for (const s of maximal) visit([...s].sort((x, y) => x - y));
  return [...seen.values()].map((vertices) => ({ vertices, value }));
}

/**
 * Lower-star filtration: each simplex enters when its highest vertex does.
 * `simplices` may be bare vertex arrays or { vertices } records; missing faces
 * are added automatically.
 */
export function lowerStar(simplices, vertexValues) {
  const all = closure(simplices.map((s) => (Array.isArray(s) ? s : s.vertices)));
  return all.map(({ vertices }) => {
    let v = -Infinity;
    for (const u of vertices) {
      const f = vertexValues[u];
      if (typeof f !== 'number' || Number.isNaN(f)) throw new RangeError(`vertex ${u} has no value`);
      if (f > v) v = f;
    }
    return { vertices, value: v };
  });
}

/** Lower-star filtration of a 1D signal on the path graph 0-1-2-...-(n-1). */
export function timeSeriesFiltration(values) {
  const simplices = [];
  for (let i = 0; i < values.length; i++) simplices.push([i]);
  for (let i = 0; i + 1 < values.length; i++) simplices.push([i, i + 1]);
  return lowerStar(simplices, values);
}

export function euclideanDistances(points) {
  const n = points.length;
  const D = Array.from({ length: n }, () => new Float64Array(n));
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const a = points[i];
      const b = points[j];
      if (a.length !== b.length) throw new Error('points have inconsistent dimension');
      let s = 0;
      for (let k = 0; k < a.length; k++) {
        const t = a[k] - b[k];
        s += t * t;
      }
      D[i][j] = D[j][i] = Math.sqrt(s);
    }
  }
  return D;
}

function checkDistanceMatrix(D) {
  const n = D.length;
  for (let i = 0; i < n; i++) {
    if (!D[i] || D[i].length !== n) throw new Error('distance matrix must be square');
    if (D[i][i] !== 0) throw new Error(`distance matrix diagonal entry ${i} is not 0`);
    for (let j = 0; j < i; j++) {
      if (D[i][j] !== D[j][i]) throw new Error(`distance matrix is not symmetric at (${i}, ${j})`);
      if (!(D[i][j] >= 0)) throw new Error(`negative or NaN distance at (${i}, ${j})`);
    }
  }
}

/**
 * Vietoris-Rips filtration of a finite metric space.
 *
 * A simplex enters at the length of its longest edge (the "diameter"
 * convention used by Ripser and GUDHI). Simplices are generated by
 * incremental clique expansion (Zomorodian 2010) over the neighbourhood graph
 * truncated at `threshold`, and only up to dimension maxDim + 1 because that
 * is all homology in degree maxDim needs.
 */
export function ripsFiltration(distances, { maxDim = 1, threshold = Infinity, maxSimplices = 2_000_000 } = {}) {
  checkDistanceMatrix(distances);
  if (!Number.isInteger(maxDim) || maxDim < 0) throw new RangeError('maxDim must be a non negative integer');
  const n = distances.length;
  const top = maxDim + 1;
  const out = [];
  const push = (vertices, value) => {
    if (out.length >= maxSimplices) {
      throw new RangeError(`Rips complex exceeds ${maxSimplices} simplices; lower threshold or maxDim`);
    }
    out.push({ vertices, value });
  };

  for (let v = 0; v < n; v++) push([v], 0);
  if (top < 1) return out;

  const upper = Array.from({ length: n }, () => []);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (distances[i][j] <= threshold) upper[i].push(j);
    }
  }
  const adjacent = upper.map((list) => new Set(list));

  const expand = (sigma, value, candidates) => {
    for (const v of candidates) {
      let tv = value;
      for (const u of sigma) if (distances[u][v] > tv) tv = distances[u][v];
      const tau = [...sigma, v];
      push(tau, tv);
      if (tau.length - 1 < top) {
        const next = candidates.filter((w) => w > v && adjacent[v].has(w));
        if (next.length) expand(tau, tv, next);
      }
    }
  };
  for (let u = 0; u < n; u++) expand([u], 0, upper[u]);
  return out;
}
