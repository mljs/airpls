import cuthillMckee from 'cuthill-mckee';

import Cholesky from './choleskySolver';

/**
 * Fit the baseline drift by iteratively changing weights of sum square error between the fitted baseline and original signals,
 * for further information about the parameters you can get the [paper of airPLS](https://github.com/zmzhang/airPLS/blob/master/airPLS_manuscript.pdf)
 * @param {Array} yData - original data
 * @param {object} [options={}] - options
 * @param {number} [options.maxIterations = 100] - maximal number of iterations if the method does not reach the stop criterion
 * @param {number} [options.factorCriterion = 0.001] - factor of the sum of absolute value of original data, to compute stop criterion
 * @param {Array} [options.weights = [1,1,...]] - initial weights vector, default each point has the same weight
 * @param {number} [options.lambda = 100] - factor of weights matrix in -> [I + lambda D'D]z = x
 * @returns {array} - list with baseline, corrected (original - baseline), iteration and error value.
 */
function airPLS(yData, options = {}) {
  let {
    maxIterations = 100,
    lambda = 100,
    factorCriterion = 0.001,
    weights = new Array(yData.length).fill(1),
    controlPoints = []
  } = options;

  var nbPoints = yData.length;
  var stopCriterion = factorCriterion * yData.reduce((sum, e) => Math.abs(e) + sum, 0);

  var { lowerTriangularNonZeros, permutationEncodedArray } = getDeltaMatrix(nbPoints, lambda);

  var sumNegDifferences = Number.MAX_SAFE_INTEGER;
  for (var iteration = 0; (iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion); iteration++) {
    let [leftHandSide, rightHandSide] = updateSystem(lowerTriangularNonZeros, yData, weights);

    let cho = Cholesky(leftHandSide, nbPoints, permutationEncodedArray);

    var baseline = cho(rightHandSide);

    sumNegDifferences = 0;

    let difference = yData.map(calculateError);

    let maxNegativeDiff = -1 * Number.MAX_SAFE_INTEGER;
    for (var i = 1, l = nbPoints - 1; i < l; i++) {
      let diff = difference[i];
      if (diff >= 0) {
        weights[i] = 0;
      } else {
        weights[i] = Math.exp(iteration * diff / sumNegDifferences);
        if (maxNegativeDiff < diff) maxNegativeDiff = diff;
      }
    }

    let value = Math.exp(iteration * maxNegativeDiff / sumNegDifferences);
    weights[0] = value;
    weights[l] = value;
    controlPoints.forEach((i) => (weights[i] = value));
  }

  return {
    corrected: yData.map((e, i) => e - baseline[i]),
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


function getDeltaMatrix(nbPoints, lambda) {
  var matrix = [];
  for (var i = 0, last = nbPoints - 1; i < last; i++) {
    matrix.push([i, i, lambda * 2]);
    matrix.push([i + 1, i, -1 * lambda]);
  }
  matrix[0][2] = lambda;
  matrix.push([last, last, lambda]);
  return { lowerTriangularNonZeros: matrix, permutationEncodedArray: cuthillMckee(matrix, nbPoints) };
}

function updateSystem(matrix, yData, weights) {
  let nbPoints = yData.length;
  var newMatrix = new Array(matrix.length);
  var newVector = new Float64Array(nbPoints);
  for (var i = 0, l = nbPoints - 1; i < l; i++) {
    let w = weights[i];
    let diag = i * 2;
    let next = diag + 1;
    newMatrix[diag] = matrix[diag].slice();
    newMatrix[next] = matrix[next].slice();
    if (w === 0) {
      newVector[i] = 0;
    } else {
      newVector[i] = yData[i] * w;
      newMatrix[diag][2] += w;
    }
  }
  newVector[l] = yData[l] * weights[l];
  newMatrix[l * 2] = matrix[l * 2].slice();
  newMatrix[l * 2][2] += weights[l];

  return [newMatrix, newVector];
}


export { airPLS as default };
