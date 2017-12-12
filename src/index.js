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
    var {A, P} = getDeltaMatrix(nbPoints, lambda, weights); //choice another name
    console.timeEnd('cov');
    return;
    console.time('decomposition');
    let cho = Cholesky.prepare(A, nbPoints, P);
    console.timeEnd('decomposition');
    console.timeEnd('prepro');
    return;
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

function getDeltaMatrix(nbPoints, lambda, weights) {
    var matrix = [];
    for (var i = 0, last = nbPoints - 1; i < last; i++) {
        matrix.push([i, i, lambda * 2 + weights[i]]);
        matrix.push([i, i + 1, -lambda]);
        matrix.push([i + 1, i, -lambda]);
    }
    matrix[0][2] = lambda + weights[0];
    matrix.push([last, last, lambda + weights[last]]);
    console.log(matrix);
    console.log(nbPoints);
    return {A: matrix, P: cuthillMckee(matrix, nbPoints)};
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