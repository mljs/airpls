import Cholesky from './choleskySolver';

import { updateSystem, getDeltaMatrix, getCloseIndex } from './utils';

/**
 * Fit the baseline drift by iteratively changing weights of sum square error between the fitted baseline and original signals,
 * for further information about the parameters you can get the [paper of airPLS](https://github.com/zmzhang/airPLS/blob/master/airPLS_manuscript.pdf)
 * @param {Array} x - x axis data useful when there ise
 * @param {Array} y - original data
 * @param {object} [options={}] - options
 * @param {number} [options.maxIterations = 100] - maximal number of iterations if the method does not reach the stop criterion
 * @param {number} [options.factorCriterion = 0.001] - factor of the sum of absolute value of original data, to compute stop criterion
 * @param {Array} [options.weights = [1,1,...]] - initial weights vector, default each point has the same weight
 * @param {number} [options.lambda = 100] - factor of weights matrix in -> [I + lambda D'D]z = x
 * @param {Array} [options.controlPoints = []] - Array of x axis values to force that baseline cross those points.
 * @param {Array} [options.baseLineZones = []] - Array of x axis values (as from - to), to force that baseline cross those zones.
 * @returns {array} - list with baseline, corrected (original - baseline), iteration and error value.
 */
function airPLS(x, y, options = {}) {
  let {
    maxIterations = 100,
    lambda = 100,
    factorCriterion = 0.001,
    weights = new Array(y.length).fill(1),
    controlPoints = [],
    baseLineZones = []
  } = options;

  if (controlPoints.length > 0) {
    controlPoints.forEach((e, i, arr) => (arr[i] = getCloseIndex(x, e)));
  }
  if (baseLineZones.length > 0) {
    baseLineZones.forEach((range) => {
      let indexFrom = getCloseIndex(x, range.from);
      let indexTo = getCloseIndex(x, range.to);
      if (indexFrom > indexTo) [indexFrom, indexTo] = [indexTo, indexFrom];
      for (let i = indexFrom; i < indexTo; i++) {
        controlPoints.push(i);
      }
    });
  }
  var nbPoints = y.length;
  var stopCriterion =
    factorCriterion * y.reduce((sum, e) => Math.abs(e) + sum, 0);

  var { lowerTriangularNonZeros, permutationEncodedArray } = getDeltaMatrix(
    nbPoints,
    lambda
  );

  var sumNegDifferences = Number.MAX_SAFE_INTEGER;
  for (
    var iteration = 0;
    iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion;
    iteration++
  ) {
    let [leftHandSide, rightHandSide] = updateSystem(
      lowerTriangularNonZeros,
      y,
      weights
    );

    let cho = Cholesky(leftHandSide, nbPoints, permutationEncodedArray);

    var baseline = cho(rightHandSide);

    sumNegDifferences = 0;

    let difference = y.map(calculateError);

    let maxNegativeDiff = -1 * Number.MAX_SAFE_INTEGER;
    for (var i = 1, l = nbPoints - 1; i < l; i++) {
      let diff = difference[i];
      if (diff >= 0) {
        weights[i] = 0;
      } else {
        weights[i] = Math.exp((iteration * diff) / sumNegDifferences);
        if (maxNegativeDiff < diff) maxNegativeDiff = diff;
      }
    }

    let value = Math.exp((iteration * maxNegativeDiff) / sumNegDifferences);
    weights[0] = value;
    weights[l] = value;
    controlPoints.forEach((i) => (weights[i] = value));
  }

  return {
    corrected: y.map((e, i) => e - baseline[i]),
    baseline,
    iteration,
    error: sumNegDifferences
  };

  function calculateError(e, i) {
    let diff = e - baseline[i];
    if (diff < 0) sumNegDifferences += diff;
    return diff;
  }
}

export { airPLS as default };
