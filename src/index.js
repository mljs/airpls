import Cholesky from 'cholesky-solve';
import cuthillMckee from 'cuthill-mckee';

function airPLS(yData, options = {}) {
    let {
        maxIterations = 100,
        lambda = 100,
        factorCriterion = 0.001,
        weights
    } = options;

    var nbPoints = yData.length;
    var stopCriterion = factorCriterion * yData.reduce((sum, e) => Math.abs(e) + sum, 0);

    if (!weights) {
        weights = new Array(nbPoints).fill(1);
    }

    var {lowerTriangularNonZeros, permutationEncodedArray} = getDeltaMatrix(nbPoints, lambda);

    var sumNegDifferences = Number.MAX_SAFE_INTEGER;
    for (var iteration = 0; (iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion); iteration++) {
        let [leftHandSide, rightHandSide] = updateSystem(lowerTriangularNonZeros, yData, weights);

        let cho = Cholesky.prepare(leftHandSide, nbPoints, permutationEncodedArray);

        var baseline = cho(rightHandSide);

        sumNegDifferences = 0;
        let difference = yData.map((e, i) => {
            let diff = e - baseline[i];
            if (diff < 0) sumNegDifferences += diff;
            return diff;
        });

        let maxNegativeDiff = -1 * Number.MAX_SAFE_INTEGER;
        for (var i = 1, l = nbPoints - 1; i < l; i++) {
            let diff = difference[i];
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
    }

    return {
        corrected: yData.map((e, i) => e - baseline[i]),
        baseline,
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
    return {lowerTriangularNonZeros: matrix, permutationEncodedArray: cuthillMckee(matrix, nbPoints)};
}

function updateSystem(matrix, yData, weights) {
    let nbPoints = yData.length;
    var newMatrix = new Array(matrix.length);
    var newVector = new Float64Array(nbPoints);
    for (var i = 0, l = nbPoints - 1; i < l; i++) {
        let w = weights[i];
        let diag = i * 2;
        let next = diag + 1;
        newMatrix[diag] = matrix[diag].slice();
        newMatrix[next] = matrix[next].slice();
        if (w === 0) {
            newVector[i] = 0;
        } else {
            newVector[i] = yData[i] * w;
            newMatrix[diag][2] += w;
        }
    }
    newVector[l] = yData[l] * weights[l];
    newMatrix[l * 2] = matrix[l * 2].slice();
    newMatrix[l * 2][2] += weights[l];

    return [newMatrix, newVector];
}

export {airPLS as default};
