import type { NumberArray } from 'cheminfo-types';
import cuthillMckee from 'cuthill-mckee';

/**
 * Updates the system matrix with new weights and computes the weighted right-hand side vector.
 * @param matrix - The lower triangular non-zeros matrix in triplet format.
 * @param y - The original data vector.
 * @param weights - The weights vector.
 * @returns A tuple of the updated matrix and the weighted right-hand side vector.
 */
export function updateSystem(
  matrix: number[][],
  y: NumberArray,
  weights: NumberArray,
): [number[][], Float64Array] {
  const nbPoints = y.length;
  const l = nbPoints - 1;
  const newMatrix = new Array<number[]>(matrix.length);
  const newVector = new Float64Array(nbPoints);
  for (let i = 0; i < l; i++) {
    const w = weights[i];
    const diag = i * 2;
    const next = diag + 1;
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
}

/**
 * Constructs the delta matrix and applies Cuthill-McKee permutation.
 * @param nbPoints - Number of data points.
 * @param lambda - Smoothing parameter.
 * @returns An object containing the lower triangular non-zeros matrix and the permutation array.
 */
export function getDeltaMatrix(
  nbPoints: number,
  lambda: number,
): {
  lowerTriangularNonZeros: number[][];
  permutationEncodedArray: number[];
} {
  const matrix: number[][] = [];
  const last = nbPoints - 1;
  for (let i = 0; i < last; i++) {
    matrix.push([i, i, lambda * 2], [i + 1, i, -1 * lambda]);
  }
  matrix[0][2] = lambda;
  matrix.push([last, last, lambda]);
  return {
    lowerTriangularNonZeros: matrix,
    permutationEncodedArray: cuthillMckee(matrix, nbPoints),
  };
}
