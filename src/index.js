import {
  xAbsolute,
  xMean,
  xMedian,
  xNoiseSanPlot,
  xStandardDeviation,
} from 'ml-spectra-processing';

import cholesky from './choleskySolver';
import { updateSystem, getDeltaMatrix, getCloseIndex } from './utils';

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
 * @param {Array<number>} [options.controlPoints = []] - Array of x axis values to force that baseline cross those points.
 * @param {Array<number>} [options.zones = []] - Array of x axis values (as from - to), to force that baseline cross those zones.
 * @returns {{corrected: Array<number>, error: number, iteration: number, baseline: Array<number>}}
 */
export default function airPLS(x, y, options = {}) {
  let {
    maxIterations = 100,
    lambda = 10,
    tolerance = 0.001,
    weights = new Array(y.length).fill(0.8),
    controlPoints = [],
    zones = [],
  } = options;
  if (controlPoints.length > 0) {
    controlPoints.forEach((e, i, arr) => (arr[i] = getCloseIndex(x, e)));
  }
  // if (zones.length > 0) {
  //   zones.forEach((range) => {
  //     let indexFrom = getCloseIndex(x, range.from);
  //     let indexTo = getCloseIndex(x, range.to);
  //     if (indexFrom > indexTo) [indexFrom, indexTo] = [indexTo, indexFrom];
  //     for (let i = indexFrom; i < indexTo; i++) {
  //       controlPoints.push(i);
  //     }
  //   });
  // }

  for (const point of controlPoints) {
    weights[point] = 2;
  }
  let baseline, iteration;
  let nbPoints = y.length;
  let l = nbPoints - 1;
  let sumNegDifferences = Number.MAX_SAFE_INTEGER;
  let stopCriterion = tolerance * y.reduce((sum, e) => Math.abs(e) + sum, 0);

  let { lowerTriangularNonZeros, permutationEncodedArray } = getDeltaMatrix(
    nbPoints,
    lambda,
  );

  let threshold = 1;
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

    let cho = cholesky(leftHandSide, nbPoints, permutationEncodedArray);

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
    console.log(
      sumNegDifferences,
      prevNegSum,
      prevNegSum / sumNegDifferences,
      stopCriterion,
      iteration,
    );

    prevNegSum = sumNegDifferences + 0;
    let maxNegativeDiff = -1 * Number.MAX_SAFE_INTEGER;
    console.log({ threshold });
    for (let i = 1; i < l; i++) {
      let diff = difference[i];

      if (iteration < 100) {
        if (Math.abs(diff) > threshold) {
          weights[i] = 0;
        } else {
          const factor = diff > 0 ? -10 : 1;
          weights[i] = Math.exp(
            (factor * (iteration * diff)) / Math.abs(sumNegDifferences),
          );
        }
      } else if (diff > 0) {
        weights[i] = 0;
      } else {
        weights[i] = Math.exp((iteration * diff) / Math.abs(sumNegDifferences));
      }

      if (diff < 0 && maxNegativeDiff < diff) maxNegativeDiff = diff;
    }

    // let value = Math.exp(
    //   (iteration * maxNegativeDiff) / Math.abs(sumNegDifferences),
    // );
    // weights[0] = 1;
    // weights[l] = 1;
    // for (const i of controlPoints) {
    //   weights[i] = value;
    // }
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
