export {
  buildFiltration,
  compareSimplices,
  closure,
  lowerStar,
  timeSeriesFiltration,
  euclideanDistances,
  ripsFiltration,
} from './filtration.js';
export { persistence, reduceBoundary, reduceCoboundary, addColumns, toPairs } from './reduce.js';
export { bottleneck, bottleneckMatching, wasserstein, hungarian, costMatrix } from './distance.js';
export {
  bettiNumber,
  bettiNumbers,
  bettiCurve,
  linspace,
  landscape,
  persistenceImage,
  persistentEntropy,
  totalPersistence,
} from './vectorize.js';
export * as datasets from './datasets.js';

import { euclideanDistances, ripsFiltration } from './filtration.js';
import { persistence } from './reduce.js';

/**
 * One call Vietoris-Rips persistence.
 *
 *   rips(points, { maxDim: 1, threshold: 2 })
 *   rips({ distances: D }, { maxDim: 2 })
 *
 * Returns the same object as persistence(): { diagrams, numSimplices, ... }
 * where diagrams[k] is the degree k diagram.
 */
export function rips(input, options = {}) {
  const { maxDim = 1, threshold = Infinity, maxSimplices, ...rest } = options;
  const D = Array.isArray(input) ? euclideanDistances(input) : input.distances;
  if (!D) throw new TypeError('rips expects an array of points or { distances }');
  const simplices = ripsFiltration(D, { maxDim, threshold, maxSimplices });
  return persistence(simplices, { ...rest, maxDim });
}
