import type { DoubleArray, NumberArray } from 'cheminfo-types';
import {
  xAbsoluteSum,
  xFindClosestIndex,
  xMultiply,
  xNoiseSanPlot,
} from 'ml-spectra-processing';

import cholesky from './choleskySolver.ts';
import { getDeltaMatrix, updateSystem } from './utils.ts';

export interface AirPLSOptions {
  /**
   * Maximal number of iterations if the method does not reach the stop criterion.
   * @default 100
   */
  maxIterations?: number;
  /**
   * Smoothing parameter. Factor of weights matrix in -> [I + lambda D'D]z = x.
   * @default 10
   */
  lambda?: number;
  /**
   * Factor of the sum of absolute value of original data, to compute stop criterion.
   * @default 0.001
   */
  tolerance?: number;
  /** Initial weights vector, default each point has the same weight. */
  weights?: NumberArray;
  /** Array of 0|1 to force the baseline to cross those points. */
  controlPoints?: NumberArray;
  /** Array of x axis ranges (as from - to), to force the baseline to cross those zones. */
  zones?: Array<{ from: number; to: number }>;
}

export interface AirPLSResult {
  /** The baseline-corrected data. */
  corrected: Float64Array;
  /** The estimated baseline. */
  baseline: number[];
  /** The number of iterations performed. */
  iteration: number;
  /** The sum of negative differences (error). */
  error: number;
}

/**
 * Fit the baseline drift by iteratively changing weights of sum square error between the fitted baseline and original signals,
 * for further information about the parameters you can get the {@link https://github.com/zmzhang/airPLS/blob/main/airPLS_manuscript.pdf paper of airPLS}.
 * @param x - X axis data useful when control points or zones are submitted.
 * @param y - Original data.
 * @param options - Options object.
 * @returns An object containing the corrected data, baseline, iteration count, and error.
 */
export default function airPLS(
  x: DoubleArray,
  y: DoubleArray,
  options: AirPLSOptions = {},
): AirPLSResult {
  const { weights, controlPoints } = getControlPoints(x, y, options);
  const { maxIterations = 100, lambda = 10, tolerance = 0.001 } = options;
  let baseline: number[] = [];
  let iteration: number;
  let sumNegDifferences = Number.MAX_SAFE_INTEGER;
  const corrected = Float64Array.from(y);
  const stopCriterion = getStopCriterion(y, tolerance);

  const { length } = y;
  const { lowerTriangularNonZeros, permutationEncodedArray } = getDeltaMatrix(
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
    const [leftHandSide, rightHandSide] = updateSystem(
      lowerTriangularNonZeros,
      y,
      weights,
    );

    const cho = cholesky(leftHandSide, length, permutationEncodedArray);
    if (cho === null) {
      throw new Error('Cholesky decomposition failed');
    }

    baseline = cho(rightHandSide);

    sumNegDifferences = applyCorrection(y, baseline, corrected);
    if (iteration === 1) {
      const { positive } = xNoiseSanPlot(corrected);
      threshold = positive;
    } else {
      const absChange = Math.abs(prevNegSum / sumNegDifferences);
      if (absChange < 1.01 && absChange > 0.99) {
        break;
      }
    }

    prevNegSum = sumNegDifferences + 0;
    const absoluteSumNegatives = Math.abs(sumNegDifferences);

    for (let i = 1; i < l; i++) {
      const absDiff = Math.abs(corrected[i]);
      const weight = Math.exp((-iteration * absDiff) / absoluteSumNegatives);
      weights[i] =
        controlPoints[i] < 1 && absDiff > threshold ? weight / 4 : weight;
    }

    weights[0] = 1;
    weights[l] = 1;
  }

  return {
    corrected,
    baseline,
    iteration,
    error: sumNegDifferences,
  };
}

function applyCorrection(
  y: DoubleArray,
  baseline: number[],
  corrected: Float64Array,
): number {
  let sumNegDifferences = 0;
  for (let i = 0; i < y.length; i++) {
    const diff = y[i] - baseline[i];
    if (diff < 0) sumNegDifferences += diff;
    corrected[i] = diff;
  }

  return sumNegDifferences;
}

function getStopCriterion(y: DoubleArray, tolerance: number): number {
  const sum = xAbsoluteSum(y);
  return tolerance * sum;
}

function getControlPoints(
  x: DoubleArray,
  y: DoubleArray,
  options: AirPLSOptions = {},
): { weights: NumberArray; controlPoints: NumberArray } {
  const { length } = x;
  const { controlPoints = Int8Array.from({ length }).fill(0) } = options;
  const { zones = [], weights = Float64Array.from({ length }).fill(0.5) } =
    options;

  if (x.length !== y.length) {
    throw new RangeError('Y should match the length with X');
  } else if (controlPoints.length !== x.length) {
    throw new RangeError('controlPoints should match the length with X');
  } else if (weights.length !== x.length) {
    throw new RangeError('weights should match the length with X');
  }

  for (const range of zones) {
    let indexFrom = xFindClosestIndex(x, range.from);
    let indexTo = xFindClosestIndex(x, range.to);
    if (indexFrom > indexTo) [indexFrom, indexTo] = [indexTo, indexFrom];
    for (let i = indexFrom; i < indexTo; i++) {
      controlPoints[i] = 1;
    }
  }

  return {
    weights:
      'controlPoints' in options || zones.length > 0
        ? xMultiply(weights, controlPoints)
        : weights,
    controlPoints,
  };
}
