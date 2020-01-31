import { toBeDeepCloseTo } from 'jest-matcher-deep-close-to';

import airPLS from '../index';

expect.extend({ toBeDeepCloseTo });

let y = [1, 1, 1, 1, 3, 6, 3, 1, 1, 1];
let x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

describe('airPLS test', () => {
  it('Small vector to find baseline', () => {
    let result = airPLS({ x, y });
    expect(result.baseline).toBeDeepCloseTo([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 2);
    expect(result.corrected).toBeDeepCloseTo([0, 0, 0, 0, 2, 5, 2, 0, 0, 0], 2);
  });
});
