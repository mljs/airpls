import airPLS from '../index';
import {toBeDeepCloseTo} from 'jest-matcher-deep-close-to';
expect.extend({toBeDeepCloseTo});

var vector = [1,1,1,1,3,6,3,1,1,1];

describe('airPLS test', () => {
    it('Small vector to find baseline', () => {
        var result = airPLS(vector);
        expect(result.baseline).toBeDeepCloseTo([1,1,1,1,1,1,1,1,1,1], 2)
        expect(result.corrected).toBeDeepCloseTo([0,0,0,0,2,5,2,0,0,0], 2)
    });
});
