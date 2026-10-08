# persistent-homology-js

Topological data analysis from first principles in plain JavaScript. Give it a point cloud (or a distance matrix, or a function on a mesh) and it tells you which loops, voids and clusters are real structure and which are noise, then lets you compare shapes with proper metrics and turn them into feature vectors for ML.

Zero dependencies. Runs in Node 18+ and in the browser.

```
$ node bin/ph.js barcode --sample circle:60 --width 40
60 points, 36050 simplices
H0  (60 bars, entropy 4.017)
  ████████████████████████████████████████▶  [0.000, inf)
  █████                                      [0.000, 0.198)
  ████                                       [0.000, 0.180)
  ...
H1  (1 bar, entropy 0.000)
       ███████████████████████████████████   [0.203, 1.674)
```

One component that never dies (H0), and exactly one long lived loop (H1). That loop is born when neighbouring samples link up at about 0.2 and only fills in at about 1.67, close to the side of the equilateral triangle inscribed in the unit circle (sqrt 3), which is when the last hole closes. That is the whole idea of persistence: features are measured by how long they survive as you zoom out.

## What is in here

| Module | What it does |
| --- | --- |
| `src/filtration.js` | Filtered simplicial complexes with validation (closed under faces, monotone values), Vietoris-Rips by incremental clique expansion, lower-star filtrations for functions on meshes and 1D signals |
| `src/reduce.js` | Boundary matrix reduction over Z/2: standard homology algorithm, the clearing (twist) optimisation, and persistent **cohomology** on the anti-transposed matrix, the trick behind Ripser. Optional representative cycles |
| `src/distance.js` | Exact bottleneck distance (binary search over candidate values plus Hopcroft-Karp) with the optimal matching, and p-Wasserstein distance via a Hungarian solver |
| `src/vectorize.js` | Betti numbers and curves, persistence landscapes, persistence images (pixel integrated Gaussians, not point samples), persistent entropy, total persistence |
| `src/datasets.js` | Seeded samplers (circle, two circles, figure eight, sphere, torus) and minimal triangulations with known homology (S², 7 vertex torus, 6 vertex RP², Klein bottle) |
| `bin/ph.js` | CLI: ASCII barcodes, JSON diagrams, and diagram distances between CSV point clouds |
| `demo/index.html` | Interactive playground: click to place points, drag the scale and watch the Rips complex, barcode and diagram update together |

## Running it

```bash
git clone https://github.com/MohammadHossinzehi/persistent-homology-js
cd persistent-homology-js
npm test                      # 30 tests, node:test, no install step
npm run bench                 # homology vs cohomology, with and without clearing

node bin/ph.js barcode points.csv --max-dim 1
node bin/ph.js barcode --sample two-circles:80 --json
node bin/ph.js distance a.csv b.csv --dim 1 --p 2
node bin/ph.js distance --sample circle:40 --sample two-circles:40
```

CSV input is one point per line, comma or whitespace separated, with an optional header row.

For the browser demo, serve the repo root with any static server (ES modules do not load from `file://`) and open `/demo/`:

```bash
npm run demo                  # or: python -m http.server 8080
# then visit http://localhost:8080/demo/
```

## Library usage

```js
import { rips, persistence, bottleneck, wasserstein, landscape, linspace, lowerStar, datasets } from './src/index.js';

// Rips persistence of a point cloud, homology up to degree 2,
// only building simplices with diameter <= 1.5
const { diagrams } = rips(points, { maxDim: 2, threshold: 1.5 });
diagrams[1]; // [{ dim: 1, birth, death, birthSimplex, deathSimplex }, ...] longest first

// Compare two shapes
bottleneck(diagramsA[1], diagramsB[1]);
wasserstein(diagramsA[1], diagramsB[1], { p: 2 });

// Fixed length features for a classifier
const grid = linspace(0, 2, 100);
const [lambda1, lambda2] = landscape(diagrams[1], grid, { levels: 2 });

// Sublevel set persistence of a height map on a triangulated grid
const result = persistence(lowerStar(triangles, heights));

// Any filtered complex you build yourself
persistence([
  { vertices: [0], value: 0 }, { vertices: [1], value: 0 }, { vertices: [2], value: 0 },
  { vertices: [0, 1], value: 1 }, { vertices: [1, 2], value: 1 }, { vertices: [0, 2], value: 1 },
  { vertices: [0, 1, 2], value: 3 },
]); // one H1 bar [1, 3)
```

Diagram functions accept either the point objects returned by `persistence` or plain `[birth, death]` pairs. Infinite deaths are `Infinity` in JS and `null` in the CLI's JSON output.

## Design notes

**Ordering.** Simplices are sorted by value, then dimension, then lexicographically. The dimension tie break guarantees a face precedes its cofaces even when they enter together, so the boundary matrix is upper triangular and the reduction is valid. The lexicographic tie break makes every run deterministic, which matters because tests compare pairings across algorithms simplex by simplex.

**Why cohomology.** The textbook algorithm reduces the boundary matrix column by column. On a Rips complex the expensive part is H1, where you reduce every triangle column even though nearly all of them end up paired with a zero length edge. Reducing the coboundary matrix instead gives the same pairs (de Silva, Morozov, Vejdemo-Johansson 2011), and clearing then runs from low dimensions up: every edge that kills a component is never reduced, and every triangle that kills a loop is never touched as a column. From `npm run bench`:

| case | simplices | homology | homology + clearing | cohomology + clearing |
| --- | ---: | ---: | ---: | ---: |
| circle n=80, H≤1 | 85,400 | 2,827,458 adds | 2,770,588 adds | **354 adds** |
| uniform box n=60, H≤1 | 36,050 | 509,639 | 499,266 | **87** |
| sphere n=70, H≤2, r≤1.2 | 16,041 | 95,121 | 64,590 | **173** |
| torus n=200, H≤2, r≤1.4 | 37,696 | 228,318 | 160,324 | **512** |

Wall time on the circle drops from about 490 ms to 45 ms, most of which is now building the complex. Cohomology without clearing is actually slower than homology (vertex coboundaries are long columns), which is a nice demonstration that the two ideas only pay off together. Representative cycles need the homology reduction, so asking for `representatives: true` switches algorithm automatically.

**Sparse Z/2 columns.** Columns are ascending arrays of row indices, so the pivot is the last entry and column addition is a merge that drops shared entries. No matrix is ever materialised.

**Rips construction.** Simplices come from clique expansion over the neighbourhood graph truncated at `threshold`, and stop at dimension `maxDim + 1` because that is all degree `maxDim` homology needs. A configurable `maxSimplices` guard throws a clear error instead of silently eating memory when someone asks for H3 of 500 points.

**Exact diagram distances.** Both distances are optimal partial matchings in which points may also go to the diagonal. The standard reduction pads each side with diagonal slots to get a square assignment problem. The bottleneck optimum is always one of the matrix entries, so a binary search over the distinct entries with a Hopcroft-Karp perfect matching test is exact, not approximate. Wasserstein solves the same matrix with costs raised to the p with a Hungarian solver. Points at infinity are matched separately: counts must agree, and pairing sorted births is optimal for any convex cost on a line.

**Persistence images** integrate the Gaussian over each pixel with the error function rather than sampling it at the pixel centre, so the total mass is right at any resolution (there is a test for that).

## Testing

`npm test` runs 30 tests with Node's built in runner. The interesting ones check the maths, not just the code paths:

* **Known topology.** Minimal triangulations of the sphere, torus, projective plane and Klein bottle give Betti numbers (1,0,1), (1,2,1), (1,1,1), (1,2,1) over Z/2. Since Z/2 cannot tell the torus from the Klein bottle, a separate test propagates orientations across the triangles to confirm which surfaces are and are not orientable.
* **Euler characteristic.** For every complex the alternating count of simplices equals the alternating sum of Betti numbers.
* **H0 equals Kruskal.** Finite H0 deaths of a Rips filtration are exactly the minimum spanning tree edge lengths, checked against an independent Kruskal with union find.
* **Four algorithms, one answer.** Homology and cohomology, with and without clearing, must produce identical pairs simplex for simplex, and cohomology with clearing must do at least 10x less work.
* **Exhaustive search.** Bottleneck, 1-Wasserstein, 2-Wasserstein and the Hungarian solver are compared to brute force over every permutation on random small inputs.
* **Stability theorem.** Moving each point by at most δ moves the Rips diagrams by at most 2δ in bottleneck distance, tested at several δ.
* **Representatives.** The cycle returned for the circle's loop has even degree at every vertex (it really is a cycle), contains its birth edge, and passes through every sample.
* **Elder rule.** Lower-star H0 of a 1D signal pairs each local minimum with the merge that kills it, younger first.

## References

* H. Edelsbrunner, J. Harer. *Computational Topology: An Introduction.* AMS, 2010.
* A. Zomorodian. Fast construction of the Vietoris-Rips complex. *Computers and Graphics*, 2010.
* C. Chen, M. Kerber. Persistent homology computation with a twist. *EuroCG*, 2011.
* V. de Silva, D. Morozov, M. Vejdemo-Johansson. Dualities in persistent (co)homology. *Inverse Problems*, 2011.
* U. Bauer. Ripser: efficient computation of Vietoris-Rips persistence barcodes. *J. Appl. Comput. Topology*, 2021.
* P. Bubenik. Statistical topological data analysis using persistence landscapes. *JMLR*, 2015.
* H. Adams et al. Persistence images: a stable vector representation of persistent homology. *JMLR*, 2017.
* D. Cohen-Steiner, H. Edelsbrunner, J. Harer. Stability of persistence diagrams. *Discrete Comput. Geom.*, 2007.

## License

MIT
