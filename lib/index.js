'use strict';

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var Cholesky = _interopDefault(require('cholesky-solve'));
var cuthillMckee = _interopDefault(require('cuthill-mckee'));

function airPLS(data, options = {}) {
    console.time('prepro');
    let {
        maxIterations = 100,
        lambda = 100,
        factorCriterion = 0.001,
        weights
    } = options;

    var nbPoints = data.length;
    var stopCriterion = factorCriterion * data.reduce((sum, e) => Math.abs(e) + sum, 0);

    if (!weights) {
        weights = new Array(nbPoints).fill(1);
    }

    console.time('cov');
    var {A, P} = getDeltaMatrix(nbPoints, lambda); //choice another name
    console.log('original', A);
    console.timeEnd('cov');
    console.timeEnd('prepro');
    var sumNegDifferences = Number.MAX_SAFE_INTEGER;
    for (var iteration = 0; (iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion); iteration++) {
        console.time('iteration: ' + iteration);
        let leftHandSide = A.slice();
        let rightHandSide = data.map((e, i) => {
            leftHandSide[i * 2][2] += weights[i];
            return [e * weights[i]];
        });
        console.log(leftHandSide);
        console.time('decomposition');
        let cho = Cholesky.prepare(leftHandSide, nbPoints, P);
        console.timeEnd('decomposition');

        console.time('solve');
        var baseline = cho(rightHandSide);
        console.timeEnd('solve');
        sumNegDifferences = 0;
        let difference = data.map((e, i) => {
            let diff = e - baseline[i];
            if (diff < 0) sumNegDifferences += diff;
            return diff;
        });

        let maxNegativeDiff = -1 * Number.MAX_SAFE_INTEGER;
        for (var i = 1, l = nbPoints - 1; i < l; i++) {
            let diff = difference[i];
            console.log(diff);
            if (diff >= 0) {
                weights[i] = 0;
            } else {
                weights[i] = Math.exp(iteration * diff / sumNegDifferences);
                if (maxNegativeDiff < diff) maxNegativeDiff = diff;
            }
        }

        let value = Math.exp(iteration * maxNegativeDiff / sumNegDifferences);
        weights[0] = value;
        weights[l] = value;
        console.timeEnd('iteration: ' + iteration);
    }
    
    return {
        corrected: data.map((e, i) => e - baseline[i]),
        baseline: baseline,
        iteration,
        error: sumNegDifferences
    };
}

function getDeltaMatrix(nbPoints, lambda) {
    var matrix = [];
    for (var i = 0, last = nbPoints - 1; i < last; i++) {
        matrix.push([i, i, lambda * 2]);
        matrix.push([i + 1, i, -1 * lambda]);
    }
    matrix[0][2] = lambda;
    matrix.push([last, last, lambda]);
    return {A: matrix, P: cuthillMckee(matrix, nbPoints)};
}

module.exports = airPLS;
