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
