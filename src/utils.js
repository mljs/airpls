import cuthillMckee from 'cuthill-mckee';

const getClosestNumber = (array = [], goal = 0) => {
  const closest = array.reduce((prev, curr) => {
    return Math.abs(curr - goal) < Math.abs(prev - goal) ? curr : prev;
  });
  return closest;
};

const getCloseIndex = (array = [], goal = 0) => {
  const closest = getClosestNumber(array, goal);
  return array.indexOf(closest);
};

const updateSystem = (matrix, y, weights) => {
  let nbPoints = y.length;
  let l = nbPoints - 1;
  let newMatrix = new Array(matrix.length);
  let newVector = new Float64Array(nbPoints);
  for (let i = 0; i < l; i++) {
    let w = weights[i];
    let diag = i * 2;
    let next = diag + 1;
    newMatrix[diag] = matrix[diag].slice();
    newMatrix[next] = matrix[next].slice();
    if (w === 0) {
      newVector[i] = 0;
    } else {
      newVector[i] = y[i] * w;
      newMatrix[diag][2] += w;
    }
  }
  newVector[l] = y[l] * weights[l];
  newMatrix[l * 2] = matrix[l * 2].slice();
  newMatrix[l * 2][2] += weights[l];

  return [newMatrix, newVector];
};

const getDeltaMatrix = (nbPoints, lambda) => {
  let matrix = [];
  let last = nbPoints - 1;
  for (let i = 0; i < last; i++) {
    matrix.push([i, i, lambda * 2]);
    matrix.push([i + 1, i, -1 * lambda]);
  }
  matrix[0][2] = lambda;
  matrix.push([last, last, lambda]);
  return {
    lowerTriangularNonZeros: matrix,
    permutationEncodedArray: cuthillMckee(matrix, nbPoints),
  };
};

export { updateSystem, getDeltaMatrix, getCloseIndex, getClosestNumber };
