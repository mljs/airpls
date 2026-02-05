import { test, expect } from 'vitest';

import airPLS from '../index';

let y = [1, 1, 1, 1, 3, 6, 3, 1, 1, 1];
let x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

test('Small vector to find baseline', () => {
  let result = airPLS(x, y);
  const expectedBaseline = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
  const expectedCorrected = [0, 0, 0, 0, 2, 5, 2, 0, 0, 0];
  expect(result.baseline).toHaveLength(expectedBaseline.length);
  expect(result.corrected).toHaveLength(expectedCorrected.length);
  result.baseline.forEach((v, i) => {
    expect(v).toBeCloseTo(expectedBaseline[i], 2);
  });
  result.corrected.forEach((v, i) => {
    expect(v).toBeCloseTo(expectedCorrected[i], 2);
  });
});
