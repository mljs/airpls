'use strict';

const {solve, Matrix, SVD} = require('ml-matrix');

function airPls(data, options = {}) {
    let {
        maxIterations = 100,
        lambda = 100,
        factorCriterion = 0.001
    } = options;

    var nbPoints = data.length;
    var stopCriterion = factorCriterion * data.reduce((sum, e) => Math.abs(e) + sum, 0);
    var dataMatrix = Matrix.columnVector(data);

    var weights = new Array(nbPoints).fill(1);
    var identityMatrix = Matrix.eye(nbPoints, nbPoints);
    var derivativeIMatrix = diffMatrix(identityMatrix);
    var covMatrix = derivativeIMatrix.transpose().mmul(derivativeIMatrix);

    var sumNegDifferences = Number.MAX_SAFE_INTEGER;
    for (var iteration = 0; (iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion); iteration++) {
        let wMatrix = Matrix.diag(weights, nbPoints, nbPoints);
        let leftHandSide = Matrix.add(wMatrix, covMatrix.clone().mul(lambda));
        let rightHandSide = wMatrix.mmul(dataMatrix); // it should be changed for array operation.

        var baseline = solve(leftHandSide, rightHandSide);
        sumNegDifferences = 0;
        let difference = data.map((e, i) => {
            let diff = e - baseline[i][0];
            if (diff < 0) sumNegDifferences += diff;
            return diff;
        });

        difference.forEach((diff, i) => {
            if (diff >= 0) {
                weights[i] = 0;
            } else {
                weights[i] = Math.exp(iteration * diff / sumNegDifferences);
            }
        });
    }
    return {
        corrected: data.map((e, i) => e - baseline[i][0]),
        baseline: baseline.to1DArray(),
        iteration
    };
}

function diffMatrix(matrix) {
    matrix = Matrix.checkMatrix(matrix);
    var A = matrix.subMatrixView(0, matrix.rows - 2, 0, matrix.columns - 1);
    var B = matrix.subMatrixView(1, matrix.rows - 1, 0, matrix.columns - 1);
    return Matrix.sub(B, A);
}

module.exports = airPls;
