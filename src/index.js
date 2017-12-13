import Cholesky from 'cholesky-solve';
import cuthillMckee from 'cuthill-mckee';

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
    console.timeEnd('cov');
    console.timeEnd('prepro');
    var sumNegDifferences = Number.MAX_SAFE_INTEGER;
    for (var iteration = 0; (iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion); iteration++) {
        console.time('iteration: ' + iteration);
        
        let [leftHandSide, rightHandSide] = updateSystem(A, data, weights);

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
    return {A: matrix, P: cuthillMckee(matrix, nbPoints)};
}

function updateSystem(matrix, vector, weights) {
    let nbPoints = vector.length;
    var newMatrix = new Array(matrix.length);
    var newVector = new Float64Array(nbPoints);
    for (var i = 0, l = nbPoints - 1; i < l; i++) {
        let w = weights[i];
        let diag = i * 2
        let next = diag + 1;
        newMatrix[diag] = matrix[diag].slice();
        newMatrix[next] = matrix[next].slice();
        if (w === 0) {
            newVector[i] = 0;
        } else {
            newVector[i] = vector[i] + w;
            newMatrix[diag][2] += w;
        }
    }
    newVector[l] = vector[l] * weights[l];
    newMatrix[l * 2] = matrix[l * 2].slice(); 
    newMatrix[l * 2][2] += weights[l];

    return [newMatrix, newVector];
}
function cholUpdate(matrix, x) {
    let r, s, c;
    let n = x.length;
    let L = matrix.clone();
    for (let i = 0; i < n; i++) {
        if (x[i] === 0) continue;
        let li = L[i][i];
        let xi = x[i]
        r = Math.sqrt(li * li + xi * xi);
        c = r / li;
        s = xi / li;
        L[i][i] = r;
        for (let j = i; j < n; i++) {
            L[j][i] = (L[j][i] + s * x[j]) / c;
            x[j] = c * x[j] - s * L[j][i];
        }
    }
    return L;
}

export {airPLS as default};