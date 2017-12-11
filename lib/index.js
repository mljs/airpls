'use strict';

var mlSparseMatrix = require('ml-sparse-matrix');

// import {CHO} from 'ml-matrix';
function airPLS(data, options = {}) {
    console.time('prepro');
    let {
        maxIterations = 100,
        lambda = 100,
        factorCriterion = 0.001
    } = options;

    var nbPoints = data.length;
    var stopCriterion = factorCriterion * data.reduce((sum, e) => Math.abs(e) + sum, 0);

    console.time('diff');
    var deltaMatrix = getDeltaMatrix(nbPoints);
    console.timeEnd('diff');
    // var wMatrix = SparseMatrix.identity(nbPoints, nbPoints);
    console.time('cov');
    var covMatrix = deltaMatrix.transpose().mmul(deltaMatrix).mul(lambda);
    console.timeEnd('cov');
    var sumNegDifferences = Number.MAX_SAFE_INTEGER;
    console.timeEnd('prepro');
    return
    for (var iteration = 0; (iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion); iteration++) {
        console.time('iteration: ' + iteration);
        let rightHandSide = covMatrix.clone();
        let vector = data.map((e, i) => {
            rightHandSide[i][i] += wMatrix[i][i];
            return [e * wMatrix[i][i]]
        });
        
        let cho = new CHO(rightHandSide);
        var baseline = cho.solve(vector);

        sumNegDifferences = 0;
        let difference = data.map((e, i) => {
            let diff = e - baseline[i][0];
            if (diff < 0) sumNegDifferences += diff;
            return diff;
        });

        let maxNegativeDiff = -1 * Number.MAX_SAFE_INTEGER;
        for (var i = 1, l = nbPoints - 1; i < l; i++) {
            let diff = difference[i];
            if (diff >= 0) {
                wMatrix.set(i, i, 0);
            } else {
                wMatrix.set(i, i, Math.exp(iteration * diff / sumNegDifferences));
                if (maxNegativeDiff < diff) maxNegativeDiff = diff;
            }
        }

        let value = Math.exp(iteration * maxNegativeDiff / sumNegDifferences);
        wMatrix.set(0, 0, value);
        wMatrix.set(l, l, value);
        console.timeEnd('iteration: ' + iteration);
    }
    
    return {
        corrected: data.map((e, i) => e - baseline[i][0]),
        baseline: baseline.to1DArray(),
        iteration
    };
}

function getDeltaMatrix(nbPoints) {
    let matrix = mlSparseMatrix.SparseMatrix.eye(nbPoints - 1, nbPoints);
    for (let i = 0, l = nbPoints-1; i < l; i++) {
        matrix.set(i, i + 1, -1);
    }
    return matrix;
}

module.exports = airPLS;
