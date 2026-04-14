import type { NumberArray } from 'cheminfo-types';
import { xFindClosestIndex, xMedian } from 'ml-spectra-processing';
import { generateSpectrum } from 'spectrum-generator';
import { expect, test } from 'vitest';

import airPLS from '../index.ts';

const y = [1, 1, 1, 1, 3, 6, 3, 1, 1, 1];
const x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

test('Small vector to find baseline', () => {
  const result = airPLS(x, y);
  const expectedBaseline = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
  const expectedCorrected = [0, 0, 0, 0, 2, 5, 2, 0, 0, 0];

  expect(result.baseline).toHaveLength(expectedBaseline.length);
  expect(result.corrected).toHaveLength(expectedCorrected.length);

  for (const [i, v] of result.baseline.entries()) {
    expect(v).toBeCloseTo(expectedBaseline[i], 2);
  }
  for (const [i, v] of result.corrected.entries()) {
    expect(v).toBeCloseTo(expectedCorrected[i], 2);
  }
});

test('cosine baseline correction with negative peak', () => {
  const from = 0;
  const to = 10;
  const nbPoints = 10240 / 2;

  const xyPeaks = [
    { x: 2, y: 2 },
    { x: 5, y: -2 },
    { x: 8, y: 2 },
  ];

  const spectrum = generateSpectrum(xyPeaks, {
    generator: {
      from,
      to,
      nbPoints,
    },
    noise: {
      distribution: 'normal',
    },
    peakOptions: {
      width: 20 / 400,
    },
  });

  const withBaseline = addCosineShape(spectrum.y, 1, 0.8);
  const xAnchorsValues = [
    0.43, 0.86, 2.875, 3.25, 3.625, 6.15, 6.5, 6.85, 9.2, 9.6,
  ];

  const controlPoints = new Int8Array(spectrum.x.length);
  for (const xAnchorValue of xAnchorsValues) {
    const index = xFindClosestIndex(spectrum.x, xAnchorValue);
    controlPoints[index] = 1;
  }

  const result = airPLS(spectrum.x, withBaseline, {
    controlPoints,
    maxIterations: 10,
    autoDownsample: false,
    lambda: 50,
    tolerance: 0.01,
  });

  const medianResult = xMedian(result.corrected);
  const medianWithBaseline = xMedian(withBaseline);

  expect(medianResult).toBeCloseTo(0, 1);
  expect(medianWithBaseline).toBeLessThan(medianResult);
});

test('cosine baseline correction with negative peak and downsampling', () => {
  const from = 0;
  const to = 10;
  const nbPoints = 10240 / 2;

  const xyPeaks = [
    { x: 2, y: 2 },
    { x: 5, y: -2 },
    { x: 8, y: 2 },
  ];

  const spectrum = generateSpectrum(xyPeaks, {
    generator: {
      from,
      to,
      nbPoints,
    },
    noise: {
      distribution: 'normal',
    },
    peakOptions: {
      width: 20 / 400,
    },
  });

  const withBaseline = addCosineShape(spectrum.y, 1, 0.8);
  const xAnchorsValues = [
    0.43, 0.86, 2.875, 3.25, 3.625, 6.15, 6.5, 6.85, 9.2, 9.6,
  ];

  const controlPoints = new Int8Array(spectrum.x.length);
  for (const xAnchorValue of xAnchorsValues) {
    const index = xFindClosestIndex(spectrum.x, xAnchorValue);
    controlPoints[index] = 1;
  }

  const result = airPLS(spectrum.x, withBaseline, {
    controlPoints,
    maxIterations: 10,
    autoDownsample: true,
    lambda: 50,
    tolerance: 0.01,
  });

  const medianResult = xMedian(result.corrected);
  const medianWithBaseline = xMedian(withBaseline);

  expect(medianResult).toBeCloseTo(0, 1);
  expect(medianWithBaseline).toBeLessThan(medianResult);
});

function addCosineShape(
  data: NumberArray,
  amplitude: number,
  frequency: number,
  phase = 0,
): number[] {
  const n = data.length;
  const result = new Array<number>(n);

  for (let i = 0; i < n; i++) {
    const cosineValue =
      amplitude * Math.cos((2 * Math.PI * frequency * i) / n + phase);
    result[i] = data[i] + cosineValue;
  }

  return result;
}
