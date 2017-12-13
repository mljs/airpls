import airPLS from '../index';
import {vector} from '../../data/vector';

describe('airPLS test', () => {
    it('Something to test', () => {
        console.time('process');
        var result = airPLS(vector, {});
        console.timeEnd('process');
        console.log(result.iteration);
    });
});
