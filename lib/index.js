'use strict';

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var cuthillMckee = _interopDefault(require('cuthill-mckee'));

function ldlSymbolic(
  n /* A and L are n-by-n, where n >= 0 */,
  Ap /* input of size n + 1, not modified */,
  Ai /* input of size nz=Ap[n], not modified */,
  Lp /* output of size n + 1, not defined on input */,
  Parent /* output of size n, not defined on input */,
  Lnz /* output of size n, not defined on input */,
  Flag /* workspace of size n, not defn. on input or output */
) {
  var i, k, p, kk, p2;

  for (k = 0; k < n; k++) {
    /* L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k) */
    Parent[k] = -1; /* parent of k is not yet known */
    Flag[k] = k; /* mark node k as visited */
    Lnz[k] = 0; /* count of nonzeros in column k of L */
    kk = k; /* kth original, or permuted, column */
    p2 = Ap[kk + 1];
    for (p = Ap[kk]; p < p2; p++) {
      /* A (i,k) is nonzero (original or permuted A) */
      i = Ai[p];

      if (i < k) {
        /* follow path from i to root of etree, stop at flagged node */
        for (; Flag[i] !== k; i = Parent[i]) {
          /* find parent of i if not yet determined */
          if (Parent[i] === -1) Parent[i] = k;
          Lnz[i]++; /* L (k,i) is nonzero */
          Flag[i] = k; /* mark i as visited */
        }
      }
    }
  }
  /* construct Lp index array from Lnz column counts */
  Lp[0] = 0;
  for (k = 0; k < n; k++) {
    Lp[k + 1] = Lp[k] + Lnz[k];
  }
}

function ldlNumeric(
  n /* A and L are n-by-n, where n >= 0 */,
  Ap /* input of size n+1, not modified */,
  Ai /* input of size nz=Ap[n], not modified */,
  Ax /* input of size nz=Ap[n], not modified */,
  Lp /* input of size n+1, not modified */,
  Parent /* input of size n, not modified */,
  Lnz /* output of size n, not defn. on input */,
  Li /* output of size lnz=Lp[n], not defined on input */,
  Lx /* output of size lnz=Lp[n], not defined on input */,
  D /* output of size n, not defined on input */,
  Y /* workspace of size n, not defn. on input or output */,
  Pattern /* workspace of size n, not defn. on input or output */,
  Flag /* workspace of size n, not defn. on input or output */
) {
  var yi, lKi;
  var i, k, p, kk, p2, len, top;
  for (k = 0; k < n; k++) {
    /* compute nonzero Pattern of kth row of L, in topological order */
    Y[k] = 0.0; /* Y(0:k) is now all zero */
    top = n; /* stack for pattern is empty */
    Flag[k] = k; /* mark node k as visited */
    Lnz[k] = 0; /* count of nonzeros in column k of L */
    kk = k; /* kth original, or permuted, column */
    p2 = Ap[kk + 1];
    for (p = Ap[kk]; p < p2; p++) {
      i = Ai[p]; /* get A(i,k) */
      if (i <= k) {
        Y[i] += Ax[p]; /* scatter A(i,k) into Y (sum duplicates) */
        for (len = 0; Flag[i] !== k; i = Parent[i]) {
          Pattern[len++] = i; /* L(k,i) is nonzero */
          Flag[i] = k; /* mark i as visited */
        }
        while (len > 0) Pattern[--top] = Pattern[--len];
      }
    }
    /* compute numerical values kth row of L (a sparse triangular solve) */
    D[k] = Y[k]; /* get D(k,k) and clear Y(k) */
    Y[k] = 0.0;
    for (; top < n; top++) {
      i = Pattern[top]; /* Pattern[top:n-1] is pattern of L(:,k) */
      yi = Y[i]; /* get and clear Y(i) */
      Y[i] = 0.0;
      p2 = Lp[i] + Lnz[i];
      for (p = Lp[i]; p < p2; p++) {
        Y[Li[p]] -= Lx[p] * yi;
      }
      lKi = yi / D[i]; /* the nonzero entry L(k,i) */
      D[k] -= lKi * yi;
      Li[p] = k; /* store L(k,i) in column form of L */
      Lx[p] = lKi;
      Lnz[i]++; /* increment count of nonzeros in col i */
    }

    if (D[k] === 0.0) return k; /* failure, D(k,k) is zero */
  }

  return n; /* success, diagonal of D is all nonzero */
}

function ldlLsolve(
  n /* L is n-by-n, where n >= 0 */,
  X /* size n. right-hand-side on input, soln. on output */,
  Lp /* input of size n+1, not modified */,
  Li /* input of size lnz=Lp[n], not modified */,
  Lx /* input of size lnz=Lp[n], not modified */
) {
  var j, p, p2;
  for (j = 0; j < n; j++) {
    p2 = Lp[j + 1];
    for (p = Lp[j]; p < p2; p++) {
      X[Li[p]] -= Lx[p] * X[j];
    }
  }
}
function ldlDsolve(
  n /* D is n-by-n, where n >= 0 */,
  X /* size n. right-hand-side on input, soln. on output */,
  D /* input of size n, not modified */
) {
  var j;
  for (j = 0; j < n; j++) {
    X[j] /= D[j];
  }
}
function ldlLTsolve(
  n /* L is n-by-n, where n >= 0 */,
  X /* size n. right-hand-side on input, soln. on output */,
  Lp /* input of size n+1, not modified */,
  Li /* input of size lnz=Lp[n], not modified */,
  Lx /* input of size lnz=Lp[n], not modified */
) {
  var j, p, p2;
  for (j = n - 1; j >= 0; j--) {
    p2 = Lp[j + 1];
    for (p = Lp[j]; p < p2; p++) {
      X[j] -= Lx[p] * X[Li[p]];
    }
  }
}

function ldlPerm(
  n /* size of X, B, and P */,
  X /* output of size n. */,
  B /* input of size n. */,
  P /* input permutation array of size n. */
) {
  var j;
  for (j = 0; j < n; j++) {
    X[j] = B[P[j]];
  }
}

function ldlPermt(
  n /* size of X, B, and P */,
  X /* output of size n. */,
  B /* input of size n. */,
  P /* input permutation array of size n. */
) {
  var j;
  for (j = 0; j < n; j++) {
    X[P[j]] = B[j];
  }
}

function prepare(M, n, P) {
  // if a permutation was specified, apply it.
  if (P) {
    var Pinv = new Array(n);

    for (let k = 0; k < n; k++) {
      Pinv[P[k]] = k;
    }

    var Mt = []; // scratch memory
    // Apply permutation. We make M into P*M*P^T
    for (let a = 0; a < M.length; ++a) {
      var ar = Pinv[M[a][0]];
      var ac = Pinv[M[a][1]];

      // we only store the upper-diagonal elements(since we assume matrix is symmetric, we only need to store these)
      // if permuted element is below diagonal, we simply transpose it.
      if (ac < ar) {
        var t = ac;
        ac = ar;
        ar = t;
      }

      Mt[a] = [];
      Mt[a][0] = ar;
      Mt[a][1] = ac;
      Mt[a][2] = M[a][2];
    }

    M = Mt; // copy scratch memory.
  } else {
    // if P argument is null, we just use an identity permutation.
    P = [];
    for (let i = 0; i < n; ++i) {
      P[i] = i;
    }
  }

  // The sparse matrix we are decomposing is A.
  // Now we shall create A from M.
  var Ap = new Array(n + 1);
  var Ai = new Array(M.length);
  var Ax = new Array(M.length);

  // count number of non-zero elements in columns.
  var LNZ = [];
  for (let i = 0; i < n; ++i) {
    LNZ[i] = 0;
  }
  for (let a = 0; a < M.length; ++a) {
    LNZ[M[a][1]]++;
  }

  Ap[0] = 0;
  for (let i = 0; i < n; ++i) {
    Ap[i + 1] = Ap[i] + LNZ[i];
  }

  var coloffset = [];
  for (let a = 0; a < n; ++a) {
    coloffset[a] = 0;
  }

  // go through all elements in M, and add them to sparse matrix A.
  for (let i = 0; i < M.length; ++i) {
    var e = M[i];
    var col = e[1];

    var adr = Ap[col] + coloffset[col];
    Ai[adr] = e[0];
    Ax[adr] = e[2];

    coloffset[col]++;
  }

  var D = new Array(n);
  var Y = new Array(n);
  var Lp = new Array(n + 1);
  var Parent = new Array(n);
  var Lnz = new Array(n);
  var Flag = new Array(n);
  var Pattern = new Array(n);
  var bp1 = new Array(n);
  var x = new Array(n);
  var d;

  ldlSymbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag);

  var Lx = new Array(Lp[n]);
  var Li = new Array(Lp[n]);

  d = ldlNumeric(n, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Pattern, Flag);

  if (d === n) {
    return function (b) {
      ldlPerm(n, bp1, b, P);
      ldlLsolve(n, bp1, Lp, Li, Lx);
      ldlDsolve(n, bp1, D);
      ldlLTsolve(n, bp1, Lp, Li, Lx);
      ldlPermt(n, x, bp1, P);

      return x;
    };
  } else {
    return null;
  }
}

const getClosestNumber = (array = [], goal = 0) => {
  const closest = array.reduce((prev, curr) => {
    return Math.abs(curr - goal) < Math.abs(prev - goal) ? curr : prev;
  });
  return closest;
};

const getCloseIndex = (array = [], goal = 0) => {
  const closest = getClosestNumber(array, goal);
  return array.indexOf(closest);
};

const updateSystem = (matrix, y, weights) => {
  let nbPoints = y.length;
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
      newVector[i] = y[i] * w;
      newMatrix[diag][2] += w;
    }
  }
  newVector[l] = y[l] * weights[l];
  newMatrix[l * 2] = matrix[l * 2].slice();
  newMatrix[l * 2][2] += weights[l];

  return [newMatrix, newVector];
};

const getDeltaMatrix = (nbPoints, lambda) => {
  var matrix = [];
  for (var i = 0, last = nbPoints - 1; i < last; i++) {
    matrix.push([i, i, lambda * 2]);
    matrix.push([i + 1, i, -1 * lambda]);
  }
  matrix[0][2] = lambda;
  matrix.push([last, last, lambda]);
  return {
    lowerTriangularNonZeros: matrix,
    permutationEncodedArray: cuthillMckee(matrix, nbPoints)
  };
};

/**
 * Fit the baseline drift by iteratively changing weights of sum square error between the fitted baseline and original signals,
 * for further information about the parameters you can get the [paper of airPLS](https://github.com/zmzhang/airPLS/blob/master/airPLS_manuscript.pdf)
 * @param {Array} y - original data
 * @param {object} [options={}] - options
 * @param {number} [options.maxIterations = 100] - maximal number of iterations if the method does not reach the stop criterion
 * @param {number} [options.factorCriterion = 0.001] - factor of the sum of absolute value of original data, to compute stop criterion
 * @param {Array} [options.weights = [1,1,...]] - initial weights vector, default each point has the same weight
 * @param {number} [options.lambda = 100] - factor of weights matrix in -> [I + lambda D'D]z = x
 * @param {Array} [options.controlPoints = []] - Array of x axis values to force that baseline cross those points.
 * @param {Array} [options.baseLineZones = []] - Array of x axis values (as from - to), to force that baseline cross those zones.
 * @returns {array} - list with baseline, corrected (original - baseline), iteration and error value.
 */
function airPLS(data, options = {}) {
  let { x, y } = data;
  let {
    maxIterations = 100,
    lambda = 100,
    factorCriterion = 0.001,
    weights = new Array(y.length).fill(1),
    controlPoints = [],
    baseLineZones = []
  } = options;

  if (controlPoints.length > 0) {
    controlPoints.forEach((e, i, arr) => (arr[i] = getCloseIndex(x, e)));
  }
  if (baseLineZones.length > 0) {
    baseLineZones.forEach((range) => {
      let indexFrom = getCloseIndex(x, range.from);
      let indexTo = getCloseIndex(x, range.to);
      if (indexFrom > indexTo) [indexFrom, indexTo] = [indexTo, indexFrom];
      for (let i = indexFrom; i < indexTo; i++) {
        controlPoints.push(i);
      }
    });
  }
  var nbPoints = y.length;
  var stopCriterion =
    factorCriterion * y.reduce((sum, e) => Math.abs(e) + sum, 0);

  var { lowerTriangularNonZeros, permutationEncodedArray } = getDeltaMatrix(
    nbPoints,
    lambda
  );

  var sumNegDifferences = Number.MAX_SAFE_INTEGER;
  for (
    var iteration = 0;
    iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion;
    iteration++
  ) {
    let [leftHandSide, rightHandSide] = updateSystem(
      lowerTriangularNonZeros,
      y,
      weights
    );

    let cho = prepare(leftHandSide, nbPoints, permutationEncodedArray);

    var baseline = cho(rightHandSide);

    sumNegDifferences = 0;

    let difference = y.map(calculateError);

    let maxNegativeDiff = -1 * Number.MAX_SAFE_INTEGER;
    for (var i = 1, l = nbPoints - 1; i < l; i++) {
      let diff = difference[i];
      if (diff >= 0) {
        weights[i] = 0;
      } else {
        weights[i] = Math.exp((iteration * diff) / sumNegDifferences);
        if (maxNegativeDiff < diff) maxNegativeDiff = diff;
      }
    }

    let value = Math.exp((iteration * maxNegativeDiff) / sumNegDifferences);
    weights[0] = value;
    weights[l] = value;
    controlPoints.forEach((i) => (weights[i] = value));
  }

  return {
    corrected: y.map((e, i) => e - baseline[i]),
    baseline,
    iteration,
    error: sumNegDifferences
  };

  function calculateError(e, i) {
    let diff = e - baseline[i];
    if (diff < 0) sumNegDifferences += diff;
    return diff;
  }
}

module.exports = airPLS;
