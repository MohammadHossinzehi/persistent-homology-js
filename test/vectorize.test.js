import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import {
  bettiNumber,
  bettiCurve,
  landscape,
  linspace,
  persistenceImage,
  persistentEntropy,
  totalPersistence,
  rips,
  datasets,
} from '../src/index.js';

const close = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) <= eps, `${a} != ${b}`);

test('Betti number counts half open intervals [birth, death)', () => {
  const dgm = [[0, 1], [0.5, 2], [1, Infinity]];
  assert.equal(bettiNumber(dgm, 0), 1);
  assert.equal(bettiNumber(dgm, 0.75), 2);
  assert.equal(bettiNumber(dgm, 1), 2);
  assert.equal(bettiNumber(dgm, 5), 1);
  assert.deepEqual(bettiCurve(dgm, [0, 0.75, 1, 5]), [1, 2, 2, 1]);
});

test('landscape of one bar is a tent of height half its length', () => {
  const grid = linspace(0, 4, 41);
  const [l1, l2] = landscape([[1, 3]], grid, { levels: 2 });
  close(Math.max(...l1), 1);
  close(l1[20], 1); // t = 2
  close(l1[10], 0); // t = 1
  close(l1[15], 0.5);
  assert.ok(l2.every((x) => x === 0));
});

test('landscape levels are ordered and nested bars stack', () => {
  const grid = linspace(0, 6, 61);
  const L = landscape([[0, 6], [1, 5], [2, 4]], grid, { levels: 3 });
  for (let g = 0; g < grid.length; g++) {
    assert.ok(L[0][g] >= L[1][g] && L[1][g] >= L[2][g]);
  }
  close(L[0][30], 3);
  close(L[1][30], 2);
  close(L[2][30], 1);
});

test('persistence image mass equals total weight when the window covers the Gaussians', () => {
  const dgm = [[0.5, 1.5], [1, 3]];
  const img = persistenceImage(dgm, { resolution: 60, sigma: 0.05, birthRange: [-1, 3], persistenceRange: [-1, 4] });
  const mass = img.flat().reduce((a, b) => a + b, 0);
  // weights are persistence / max persistence = 0.5 and 1
  close(mass, 1.5, 1e-5);
  // the brighter pixel sits near (birth 1, persistence 2)
  let best = [0, 0, -1];
  img.forEach((row, r) => row.forEach((v, c) => { if (v > best[2]) best = [r, c, v]; }));
  const pers = -1 + ((best[0] + 0.5) * 5) / 60;
  const birth = -1 + ((best[1] + 0.5) * 4) / 60;
  close(pers, 2, 0.1);
  close(birth, 1, 0.1);
});

test('persistent entropy and total persistence', () => {
  close(persistentEntropy([[0, 1], [0, 1]]), Math.log(2));
  close(persistentEntropy([[0, 1]]), 0);
  close(persistentEntropy([[0, 1], [5, Infinity]]), 0);
  close(totalPersistence([[0, 1], [1, 4], [0, Infinity]], 2), 10);
  // a clean circle concentrates persistence in one bar, so H1 entropy is low
  const ring = rips(datasets.circle(40, { noise: 0.02 }), { maxDim: 1 }).diagrams[1];
  const blob = rips(datasets.uniformBox(40, 2, { seed: 3 }), { maxDim: 1 }).diagrams[1];
  assert.ok(persistentEntropy(ring) < persistentEntropy(blob));
});

test('CLI prints a barcode and a JSON diagram', () => {
  const cli = fileURLToPath(new URL('../bin/ph.js', import.meta.url));
  const text = execFileSync(process.execPath, [cli, 'barcode', '--sample', 'circle:40', '--max-dim', '1'], { encoding: 'utf8' });
  assert.match(text, /H0/);
  assert.match(text, /H1/);
  const json = JSON.parse(
    execFileSync(process.execPath, [cli, 'barcode', '--sample', 'circle:40', '--json'], { encoding: 'utf8' }),
  );
  assert.equal(json.diagrams[1].length >= 1, true);
  assert.equal(json.diagrams[0][0][1], null); // Infinity is serialised as null
  const dist = execFileSync(process.execPath, [cli, 'distance', '--sample', 'circle:30', '--sample', 'two-circles:30'], {
    encoding: 'utf8',
  });
  assert.match(dist, /bottleneck/);
});
