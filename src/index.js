import { xMultiply, xNoiseSanPlot } from 'ml-spectra-processing';

import cholesky from './choleskySolver';
import { updateSystem, getDeltaMatrix, getCloseIndex } from './utils';

function getControlPoints(x, y, options = {}) {
  const { length } = x;
  let { controlPoints = Int8Array.from({ length }).fill(0) } = options;
  const { zones = [], weights = Float64Array.from({ length }).fill(1) } =
    options;

  if (x.length !== y.length) {
    throw new RangeError('Y should match the length with X');
  } else if (controlPoints.length !== x.length) {
    throw new RangeError('controlPoints should match the length with X');
  } else if (weights.length !== x.length) {
    throw new RangeError('weights should match the length with X');
  }

  zones.forEach((range) => {
    let indexFrom = getCloseIndex(x, range.from);
    let indexTo = getCloseIndex(x, range.to);
    if (indexFrom > indexTo) [indexFrom, indexTo] = [indexTo, indexFrom];
    for (let i = indexFrom; i < indexTo; i++) {
      controlPoints[i] = 1;
    }
  });

  return {
    weights:
      'controlPoints' in options || zones.length > 0
        ? xMultiply(weights, controlPoints)
        : weights,
    controlPoints,
  };
}

/**
 * Fit the baseline drift by iteratively changing weights of sum square error between the fitted baseline and original signals,
 * for further information about the parameters you can get the [paper of airPLS](https://github.com/zmzhang/airPLS/blob/master/airPLS_manuscript.pdf)
 * @param {Array<number>} x - x axis data useful when control points or zones are submitted
 * @param {Array<number>} y - Original data
 * @param {object} [options={}] - Options object
 * @param {number} [options.maxIterations = 100] - Maximal number of iterations if the method does not reach the stop criterion
 * @param {number} [options.tolerance = 0.001] - Factor of the sum of absolute value of original data, to compute stop criterion
 * @param {Array<number>} [options.weights = [1,1,...]] - Initial weights vector, default each point has the same weight
 * @param {number} [options.lambda = 100] - Factor of weights matrix in -> [I + lambda D'D]z = x
 * @param {Array<number>} [options.controlPoints = []] - Array of 0|1 to force the baseline cross those points.
 * @param {Array<number>} [options.zones = []] - Array of x axis values (as from - to), to force that baseline cross those zones.
 * @returns {{corrected: Array<number>, error: number, iteration: number, baseline: Array<number>}}
 */
export default function airPLS(x, y, options = {}) {
  const { weights, controlPoints } = getControlPoints(x, y, options);
  let { maxIterations = 100, lambda = 10, tolerance = 0.001 } = options;

  let baseline, iteration;
  let sumNegDifferences = Number.MAX_SAFE_INTEGER;
  let stopCriterion = tolerance * y.reduce((sum, e) => Math.abs(e) + sum, 0);

  const { length } = y;
  let { lowerTriangularNonZeros, permutationEncodedArray } = getDeltaMatrix(
    length,
    lambda,
  );

  let threshold = 1;
  const l = length - 1;
  let prevNegSum = Number.MAX_SAFE_INTEGER;
  for (
    iteration = 0;
    iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion;
    iteration++
  ) {
    let [leftHandSide, rightHandSide] = updateSystem(
      lowerTriangularNonZeros,
      y,
      weights,
    );

    let cho = cholesky(leftHandSide, length, permutationEncodedArray);

    baseline = cho(rightHandSide);

    sumNegDifferences = 0;
    let difference = y.map(calculateError);
    if (iteration === 1) {
      const { positive } = xNoiseSanPlot(difference);
      threshold = positive;
    } else {
      const absChange = Math.abs(prevNegSum / sumNegDifferences);
      if (absChange < 1.01 && absChange > 0.99) {
        break;
      }
    }

    prevNegSum = sumNegDifferences + 0;
    let maxNegativeDiff = -1 * Number.MAX_SAFE_INTEGER;
    for (let i = 1; i < l; i++) {
      let diff = difference[i];
      if (controlPoints[i] < 1 && Math.abs(diff) > threshold) {
        weights[i] = 0;
      } else {
        const factor = diff > 0 ? -1 : 1;
        weights[i] = Math.exp(
          (factor * (iteration * diff)) / Math.abs(sumNegDifferences),
        );
      }

      if (diff < 0 && maxNegativeDiff < diff) maxNegativeDiff = diff;
    }

    weights[0] = 1;
    weights[l] = 1;
  }

  return {
    corrected: y.map((e, i) => e - baseline[i]),
    baseline,
    iteration,
    error: sumNegDifferences,
  };

  function calculateError(e, i) {
    let diff = e - baseline[i];
    if (diff < 0) sumNegDifferences += diff;
    return diff;
  }
}
