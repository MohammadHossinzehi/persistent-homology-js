#!/usr/bin/env node
// Command line front end.
//
//   ph barcode points.csv [--max-dim 1] [--threshold 2] [--json] [--width 60]
//   ph barcode --sample circle:80
//   ph distance a.csv b.csv [--dim 1] [--p 2]
//
// CSV input: one point per line, comma or whitespace separated numbers. A
// non numeric first line is treated as a header and skipped.

import { readFileSync } from 'node:fs';
import { rips, bottleneck, wasserstein, persistentEntropy, datasets } from '../src/index.js';

const SAMPLES = {
  circle: (n) => datasets.circle(n, { noise: 0.03 }),
  'two-circles': (n) => datasets.twoCircles(n, { noise: 0.03 }),
  'figure-eight': (n) => datasets.figureEight(n, { noise: 0.01 }),
  sphere: (n) => datasets.sphere(n),
  torus: (n) => datasets.torus(n),
  box: (n) => datasets.uniformBox(n, 2),
};

function usage(code = 0) {
  const msg = `usage:
  ph barcode <points.csv | --sample NAME:N> [--max-dim K] [--threshold T] [--json] [--width W]
  ph distance <a.csv | --sample NAME:N> <b.csv | --sample NAME:N> [--dim K] [--p P] [--threshold T]

samples: ${Object.keys(SAMPLES).join(', ')}`;
  (code ? process.stderr : process.stdout).write(msg + '\n');
  process.exit(code);
}

function parseArgs(argv) {
  const opts = { inputs: [], maxDim: 1, threshold: Infinity, json: false, width: 60, dim: 1, p: 2 };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    const next = () => {
      if (i + 1 >= argv.length) usage(1);
      return argv[++i];
    };
    if (a === '--sample') opts.inputs.push({ sample: next() });
    else if (a === '--max-dim') opts.maxDim = Number(next());
    else if (a === '--threshold') opts.threshold = Number(next());
    else if (a === '--json') opts.json = true;
    else if (a === '--width') opts.width = Number(next());
    else if (a === '--dim') opts.dim = Number(next());
    else if (a === '--p') {
      const v = next();
      opts.p = v === 'inf' ? Infinity : Number(v);
    }
    else if (a === '-h' || a === '--help') usage(0);
    else if (a.startsWith('--')) {
      process.stderr.write(`unknown option ${a}\n`);
      usage(1);
    } else opts.inputs.push({ file: a });
  }
  return opts;
}

function loadPoints(input) {
  if (input.sample) {
    const [name, count = '60'] = input.sample.split(':');
    const make = SAMPLES[name];
    if (!make) throw new Error(`unknown sample "${name}"`);
    return make(Number(count));
  }
  const lines = readFileSync(input.file, 'utf8').split(/\r?\n/).filter((l) => l.trim());
  const rows = lines.map((l) => l.trim().split(/[\s,;]+/).map(Number));
  if (rows.length && rows[0].some(Number.isNaN)) rows.shift();
  rows.forEach((r, i) => {
    if (r.some(Number.isNaN)) throw new Error(`${input.file}: line ${i + 1} is not numeric`);
  });
  return rows;
}

function drawBarcode(diagrams, width) {
  const finiteMax = Math.max(
    1e-9,
    ...diagrams.flat().map((p) => (p.death === Infinity ? p.birth : p.death)),
  );
  const scale = (x) => Math.round((Math.min(x, finiteMax) / finiteMax) * width);
  const lines = [];
  diagrams.forEach((dgm, k) => {
    lines.push(`H${k}  (${dgm.length} bar${dgm.length === 1 ? '' : 's'}, entropy ${persistentEntropy(dgm).toFixed(3)})`);
    for (const p of dgm.slice(0, 12)) {
      const a = scale(p.birth);
      const b = p.death === Infinity ? width : Math.max(a + 1, scale(p.death));
      const bar = ' '.repeat(a) + '█'.repeat(b - a) + (p.death === Infinity ? '▶' : '');
      const label = `[${p.birth.toFixed(3)}, ${p.death === Infinity ? 'inf' : p.death.toFixed(3)})`;
      lines.push(`  ${bar.padEnd(width + 2)} ${label}`);
    }
    if (dgm.length > 12) lines.push(`  ... ${dgm.length - 12} shorter bars omitted`);
  });
  lines.push(`  ${'0'.padEnd(width)}${finiteMax.toFixed(3)}`);
  return lines.join('\n');
}

function main() {
  const [command, ...rest] = process.argv.slice(2);
  if (!command || command === '-h' || command === '--help') usage(0);
  const opts = parseArgs(rest);

  if (command === 'barcode') {
    if (opts.inputs.length !== 1) usage(1);
    const points = loadPoints(opts.inputs[0]);
    const result = rips(points, { maxDim: opts.maxDim, threshold: opts.threshold });
    if (opts.json) {
      const out = {
        points: points.length,
        simplices: result.numSimplices,
        diagrams: result.diagrams.map((d) => d.map((p) => [p.birth, p.death === Infinity ? null : p.death])),
      };
      process.stdout.write(JSON.stringify(out) + '\n');
    } else {
      process.stdout.write(`${points.length} points, ${result.numSimplices} simplices\n`);
      process.stdout.write(drawBarcode(result.diagrams, opts.width) + '\n');
    }
  } else if (command === 'distance') {
    if (opts.inputs.length !== 2) usage(1);
    const [a, b] = opts.inputs.map((inp) =>
      rips(loadPoints(inp), { maxDim: opts.dim, threshold: opts.threshold }).diagrams[opts.dim],
    );
    process.stdout.write(`H${opts.dim} bottleneck   ${bottleneck(a, b).toFixed(6)}\n`);
    process.stdout.write(`H${opts.dim} wasserstein-${opts.p} ${wasserstein(a, b, { p: opts.p }).toFixed(6)}\n`);
  } else {
    usage(1);
  }
}

try {
  main();
} catch (err) {
  process.stderr.write(`error: ${err.message}\n`);
  process.exit(1);
}
