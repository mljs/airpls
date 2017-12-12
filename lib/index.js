'use strict';

var mlMatrix = require('ml-matrix');
var mlSparseMatrix = require('ml-sparse-matrix');

function airPLS(data, options = {}) {
    console.time('prepro');
    let {
        maxIterations = 100,
        lambda = 100,
        factorCriterion = 0.001
    } = options;

    var nbPoints = data.length;
    var stopCriterion = factorCriterion * data.reduce((sum, e) => Math.abs(e) + sum, 0);

    var wMatrix = mlSparseMatrix.SparseMatrix.identity(nbPoints, nbPoints);
    console.time('cov');
    var covMatrix = getDeltaMatrix(nbPoints, lambda);
    console.timeEnd('cov');
    return;
    console.time('decomposition');
    let cho = new mlMatrix.CHO(covMatrix);
    console.timeEnd('decomposition');
    console.timeEnd('prepro');
    console.time('decomposition');
    let cho = new mlMatrix.CHO(covMatrix);
    console.timeEnd('decomposition');
    console.timeEnd('prepro');
    console.time('decomposition');
    let cho = new mlMatrix.CHO(covMatrix);
    console.timeEnd('decomposition');
    console.timeEnd('prepro');
    console.time('decomposition');
    let cho = new mlMatrix.CHO(covMatrix);
    console.timeEnd('decomposition');
    console.timeEnd('prepro');
    
    var sumNegDifferences = Number.MAX_SAFE_INTEGER;
    for (var iteration = 0; (iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion); iteration++) {
        console.time('iteration: ' + iteration);

        let vector = data.map((e, i) => {
            rightHandSide[i][i] += wMatrix[i][i];
            return [e * wMatrix[i][i]]
        });
        
        
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

function getDeltaMatrix(nbPoints, lambda) {
    let matrix = new mlSparseMatrix.SparseMatrix(nbPoints, nbPoints);
    for (var i = 0; i < nbPoints; i++) {
        matrix.set(i, i, lambda * 2);
        matrix.set(i, i + 1, -lambda);
        matrix.set(i + 1, i, -lambda);
    }
    let last = nbPoints - 1;
    matrix.set(0, 0, lambda);
    matrix.set(last, last, lambda);
    return matrix;
}

module.exports = airPLS;
