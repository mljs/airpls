import type { DoubleArray, NumberArray } from 'cheminfo-types';
import {
  xAbsoluteSum,
  xFindClosestIndex,
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
  /**
   * Enable automatic downsampling for large datasets to speed up computation.
   * @default false
   */
  autoDownsample?: boolean;
  /**
   * Maximum resolution (number of points) before downsampling is applied.
   * Only used if autoDownsample is true.
   * @default 5000
   */
  maxResolution?: number;
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
  const { autoDownsample = true, maxResolution = 5000 } = options;

  // Check if downsampling should be applied
  const shouldDownsample = autoDownsample && y.length > maxResolution;
  let xWork = x;
  let yWork = y;
  let downsampleFactor = 1;
  let optionsWork = options;

  if (shouldDownsample) {
    downsampleFactor = getDownsampleFactor(y.length, maxResolution);
    yWork = averagePool(y, downsampleFactor);
    xWork = decimateIndices(x, downsampleFactor);

    // Downsample controlPoints if provided, to match downsampled x and y
    if (options.controlPoints) {
      const { controlPoints } = options;

      const downsampledControlPoints = new Int8Array(xWork.length);
      for (let i = 0; i < x.length; i++) {
        if (controlPoints[i] > 0) {
          const closestIndex = xFindClosestIndex(xWork, x[i]);
          downsampledControlPoints[closestIndex] = 1;
        }
      }

      optionsWork = {
        ...options,
        controlPoints: downsampledControlPoints,
      };
    }
  }

  const { weights, controlPoints } = getControlPoints(
    xWork,
    yWork,
    optionsWork,
  );
  const { maxIterations = 100, lambda = 10, tolerance = 0.001 } = options;
  let baseline: number[] = [];
  let iteration: number;
  let sumNegDifferences = Number.MAX_SAFE_INTEGER;
  const corrected = Float64Array.from(yWork);
  const stopCriterion = getStopCriterion(yWork, tolerance);

  const { length } = yWork;
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
      yWork,
      weights,
    );

    const cho = cholesky(leftHandSide, length, permutationEncodedArray);
    if (cho === null) {
      throw new Error('Cholesky decomposition failed');
    }

    baseline = cho(rightHandSide);

    sumNegDifferences = applyCorrection(yWork, baseline, corrected);
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

    for (let i = 1; i < l; i++) {
      const diff = corrected[i];
      if (controlPoints[i] < 1 && Math.abs(diff) > threshold) {
        weights[i] = 0;
      } else {
        const factor = diff > 0 ? -1 : 1;
        weights[i] = Math.exp(
          (factor * (iteration * diff)) / Math.abs(sumNegDifferences),
        );
      }
    }

    weights[0] = 1;
    weights[l] = 1;
  }

  // Interpolate results back to original resolution if downsampling was applied
  let finalBaseline = baseline;
  let finalCorrected = corrected;

  if (shouldDownsample) {
    finalBaseline = interpolateLinear(xWork, baseline, x);
    finalCorrected = Float64Array.from(
      y.map((val, i) => val - finalBaseline[i]),
    );
  }

  return {
    corrected: finalCorrected,
    baseline: finalBaseline,
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
  const { zones = [], weights = Float64Array.from({ length }).fill(1) } =
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
    weights,
    controlPoints,
  };
}

/**
 * Calculate the downsampling factor to reduce data to target resolution.
 * @param originalLength - Original data length.
 * @param targetResolution - Target number of points.
 * @returns Downsampling factor.
 */
function getDownsampleFactor(
  originalLength: number,
  targetResolution: number,
): number {
  return Math.max(1, Math.ceil(originalLength / targetResolution));
}

/**
 * Downsample by averaging consecutive points (average pooling).
 * @param arr - Input array.
 * @param poolSize - Number of consecutive points to average.
 * @returns Downsampled array.
 */
function averagePool(arr: DoubleArray, poolSize: number): Float64Array {
  if (poolSize <= 1) return Float64Array.from(arr);

  const result: number[] = [];
  for (let i = 0; i < arr.length; i += poolSize) {
    let sum = 0;
    const endIdx = Math.min(i + poolSize, arr.length);
    for (let j = i; j < endIdx; j++) {
      sum += arr[j];
    }
    result.push(sum / (endIdx - i));
  }
  return Float64Array.from(result);
}

/**
 * Downsample by keeping every N-th index (decimation for x-axis).
 * @param arr - Input x-axis array.
 * @param factor - Decimation factor.
 * @returns Decimated array.
 */
function decimateIndices(arr: DoubleArray, factor: number): Float64Array {
  if (factor <= 1) return Float64Array.from(arr);

  const result: number[] = [];
  for (let i = 0; i < arr.length; i += factor) {
    result.push(arr[i]);
  }
  return Float64Array.from(result);
}

/**
 * Interpolate values using linear interpolation.
 * @param xSparse - Sparse x-axis values (downsampled).
 * @param ySparse - Sparse y values (downsampled).
 * @param xTarget - Target x-axis values (original resolution).
 * @returns Interpolated y values.
 */
function interpolateLinear(
  xSparse: DoubleArray,
  ySparse: DoubleArray,
  xTarget: DoubleArray,
): number[] {
  const result = new Array(xTarget.length);

  for (let i = 0; i < xTarget.length; i++) {
    const targetX = xTarget[i];

    // Find surrounding points using binary search
    let left = 0;
    let right = xSparse.length - 1;

    while (left < right) {
      const mid = Math.floor((left + right) / 2);
      if (xSparse[mid] < targetX) {
        left = mid + 1;
      } else {
        right = mid;
      }
    }

    // Handle boundary cases
    if (left === 0) {
      if (xSparse[0] === targetX) {
        result[i] = ySparse[0];
      } else if (xSparse.length === 1) {
        result[i] = ySparse[0];
      } else {
        // Linear interpolation between first two points
        const t = (targetX - xSparse[0]) / (xSparse[1] - xSparse[0]);
        result[i] = ySparse[0] * (1 - t) + ySparse[1] * t;
      }
    } else if (left >= xSparse.length - 1) {
      result[i] = ySparse[xSparse.length - 1];
    } else {
      // Linear interpolation between two surrounding points
      const leftIdx = left - 1;
      const rightIdx = left;
      const t =
        (targetX - xSparse[leftIdx]) / (xSparse[rightIdx] - xSparse[leftIdx]);
      result[i] = ySparse[leftIdx] * (1 - t) + ySparse[rightIdx] * t;
    }
  }

  return result;
}
