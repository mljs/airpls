'use strict';

var cuthillMckee = require('cuthill-mckee');

function _interopDefaultLegacy(e) { return e && typeof e === 'object' && 'default' in e ? e : { 'default': e }; }

var cuthillMckee__default = /*#__PURE__*/_interopDefaultLegacy(cuthillMckee);

// Based on https://github.com/scijs/cholesky-solve

/*
The MIT License (MIT)

Copyright (c) 2013 Eric Arnebäck

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

function ldlSymbolic(
  n /* A and L are n-by-n, where n >= 0 */,
  Ap /* input of size n + 1, not modified */,
  Ai /* input of size nz=Ap[n], not modified */,
  Lp /* output of size n + 1, not defined on input */,
  Parent /* output of size n, not defined on input */,
  Lnz /* output of size n, not defined on input */,
  Flag /* workspace of size n, not defn. on input or output */,
) {
  let i, k, p, kk, p2;

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
  Flag /* workspace of size n, not defn. on input or output */,
) {
  let yi, lKi;
  let i, k, p, kk, p2, len, top;
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
  Lx /* input of size lnz=Lp[n], not modified */,
) {
  let j, p, p2;
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
  D /* input of size n, not modified */,
) {
  let j;
  for (j = 0; j < n; j++) {
    X[j] /= D[j];
  }
}
function ldlLTsolve(
  n /* L is n-by-n, where n >= 0 */,
  X /* size n. right-hand-side on input, soln. on output */,
  Lp /* input of size n+1, not modified */,
  Li /* input of size lnz=Lp[n], not modified */,
  Lx /* input of size lnz=Lp[n], not modified */,
) {
  let j, p, p2;
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
  P /* input permutation array of size n. */,
) {
  let j;
  for (j = 0; j < n; j++) {
    X[j] = B[P[j]];
  }
}

function ldlPermt(
  n /* size of X, B, and P */,
  X /* output of size n. */,
  B /* input of size n. */,
  P /* input permutation array of size n. */,
) {
  let j;
  for (j = 0; j < n; j++) {
    X[P[j]] = B[j];
  }
}

function prepare(M, n, P) {
  // if a permutation was specified, apply it.
  if (P) {
    let Pinv = new Array(n);

    for (let k = 0; k < n; k++) {
      Pinv[P[k]] = k;
    }

    let Mt = []; // scratch memory
    // Apply permutation. We make M into P*M*P^T
    for (let a = 0; a < M.length; ++a) {
      let ar = Pinv[M[a][0]];
      let ac = Pinv[M[a][1]];

      // we only store the upper-diagonal elements(since we assume matrix is symmetric, we only need to store these)
      // if permuted element is below diagonal, we simply transpose it.
      if (ac < ar) {
        let t = ac;
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
  let Ap = new Array(n + 1);
  let Ai = new Array(M.length);
  let Ax = new Array(M.length);

  // count number of non-zero elements in columns.
  let LNZ = [];
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

  let coloffset = [];
  for (let a = 0; a < n; ++a) {
    coloffset[a] = 0;
  }

  // go through all elements in M, and add them to sparse matrix A.
  for (let i = 0; i < M.length; ++i) {
    let e = M[i];
    let col = e[1];

    let adr = Ap[col] + coloffset[col];
    Ai[adr] = e[0];
    Ax[adr] = e[2];

    coloffset[col]++;
  }

  let D = new Array(n);
  let Y = new Array(n);
  let Lp = new Array(n + 1);
  let Parent = new Array(n);
  let Lnz = new Array(n);
  let Flag = new Array(n);
  let Pattern = new Array(n);
  let bp1 = new Array(n);
  let x = new Array(n);
  let d;

  ldlSymbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag);

  let Lx = new Array(Lp[n]);
  let Li = new Array(Lp[n]);

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
  let l = nbPoints - 1;
  let newMatrix = new Array(matrix.length);
  let newVector = new Float64Array(nbPoints);
  for (let i = 0; i < l; i++) {
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
  let matrix = [];
  let last = nbPoints - 1;
  for (let i = 0; i < last; i++) {
    matrix.push([i, i, lambda * 2]);
    matrix.push([i + 1, i, -1 * lambda]);
  }
  matrix[0][2] = lambda;
  matrix.push([last, last, lambda]);
  return {
    lowerTriangularNonZeros: matrix,
    permutationEncodedArray: cuthillMckee__default['default'](matrix, nbPoints),
  };
};

/**
 * Fit the baseline drift by iteratively changing weights of sum square error between the fitted baseline and original signals,
 * for further information about the parameters you can get the [paper of airPLS](https://github.com/zmzhang/airPLS/blob/main/airPLS_manuscript.pdf)
 * @param {Array<number>} x - x axis data useful when control points or zones are submitted
 * @param {Array<number>} y - Original data
 * @param {object} [options={}] - Options object
 * @param {number} [options.maxIterations = 100] - Maximal number of iterations if the method does not reach the stop criterion
 * @param {number} [options.factorCriterion = 0.001] - Factor of the sum of absolute value of original data, to compute stop criterion
 * @param {Array<number>} [options.weights = [1,1,...]] - Initial weights vector, default each point has the same weight
 * @param {number} [options.lambda = 100] - Factor of weights matrix in -> [I + lambda D'D]z = x
 * @param {Array<number>} [options.controlPoints = []] - Array of x axis values to force that baseline cross those points.
 * @param {Array<number>} [options.baseLineZones = []] - Array of x axis values (as from - to), to force that baseline cross those zones.
 * @returns {{corrected: Array<number>, error: number, iteration: number, baseline: Array<number>}}
 */
function airPLS(x, y, options = {}) {
  let {
    maxIterations = 100,
    lambda = 100,
    factorCriterion = 0.001,
    weights = new Array(y.length).fill(1),
    controlPoints = [],
    baseLineZones = [],
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

  let baseline, iteration;
  let nbPoints = y.length;
  let l = nbPoints - 1;
  let sumNegDifferences = Number.MAX_SAFE_INTEGER;
  let stopCriterion =
    factorCriterion * y.reduce((sum, e) => Math.abs(e) + sum, 0);

  let { lowerTriangularNonZeros, permutationEncodedArray } = getDeltaMatrix(
    nbPoints,
    lambda,
  );

  for (
    iteration = 0;
    iteration < maxIterations && Math.abs(sumNegDifferences) > stopCriterion;
    iteration++
  ) {
    let [leftHandSide, rightHandSide] = updateSystem(
      lowerTriangularNonZeros,
      y,
      weights,
    );

    let cho = prepare(leftHandSide, nbPoints, permutationEncodedArray);

    baseline = cho(rightHandSide);

    sumNegDifferences = 0;

    let difference = y.map(calculateError);

    let maxNegativeDiff = -1 * Number.MAX_SAFE_INTEGER;
    for (let i = 1; i < l; i++) {
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
    error: sumNegDifferences,
  };

  function calculateError(e, i) {
    let diff = e - baseline[i];
    if (diff < 0) sumNegDifferences += diff;
    return diff;
  }
}

module.exports = airPLS;
