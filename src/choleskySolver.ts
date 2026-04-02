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

import type { NumberArray } from 'cheminfo-types';

function ldlSymbolic(
  n: number,
  Ap: number[],
  Ai: number[],
  Lp: number[],
  Parent: number[],
  Lnz: number[],
  Flag: number[],
): void {
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
  n: number,
  Ap: number[],
  Ai: number[],
  Ax: number[],
  Lp: number[],
  Parent: number[],
  Lnz: number[],
  Li: number[],
  Lx: number[],
  D: number[],
  Y: number[],
  Pattern: number[],
  Flag: number[],
): number {
  let yi, lKi;
  let i, k, p, kk, p2, len, top;
  for (k = 0; k < n; k++) {
    /* compute nonzero Pattern of kth row of L, in topological order */
    Y[k] = 0; /* Y(0:k) is now all zero */
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
    Y[k] = 0;
    for (; top < n; top++) {
      i = Pattern[top]; /* Pattern[top:n-1] is pattern of L(:,k) */
      yi = Y[i]; /* get and clear Y(i) */
      Y[i] = 0;
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

    if (D[k] === 0) return k; /* failure, D(k,k) is zero */
  }

  return n; /* success, diagonal of D is all nonzero */
}

function ldlLsolve(
  n: number,
  X: number[],
  Lp: number[],
  Li: number[],
  Lx: number[],
): void {
  let j, p, p2;
  for (j = 0; j < n; j++) {
    p2 = Lp[j + 1];
    for (p = Lp[j]; p < p2; p++) {
      X[Li[p]] -= Lx[p] * X[j];
    }
  }
}

function ldlDsolve(n: number, X: number[], D: number[]): void {
  let j;
  for (j = 0; j < n; j++) {
    X[j] /= D[j];
  }
}

function ldlLTsolve(
  n: number,
  X: number[],
  Lp: number[],
  Li: number[],
  Lx: number[],
): void {
  let j, p, p2;
  for (j = n - 1; j >= 0; j--) {
    p2 = Lp[j + 1];
    for (p = Lp[j]; p < p2; p++) {
      X[j] -= Lx[p] * X[Li[p]];
    }
  }
}

function ldlPerm(n: number, X: number[], B: NumberArray, P: number[]): void {
  let j;
  for (j = 0; j < n; j++) {
    X[j] = B[P[j]];
  }
}

function ldlPermt(n: number, X: number[], B: number[], P: number[]): void {
  let j;
  for (j = 0; j < n; j++) {
    X[P[j]] = B[j];
  }
}

/**
 * Prepares a sparse matrix for Cholesky decomposition and returns a solver function.
 * @param M - Sparse matrix in triplet format (array of [row, col, value] entries).
 * @param n - Size of the matrix (n-by-n).
 * @param P - Optional permutation array of size n.
 * @returns A solver function that takes a right-hand side vector and returns the solution, or null if decomposition fails.
 */
function prepare(
  M: number[][],
  n: number,
  P?: number[],
): ((b: NumberArray) => number[]) | null {
  // if a permutation was specified, apply it.
  if (P) {
    const Pinv = new Array<number>(n);

    for (let k = 0; k < n; k++) {
      Pinv[P[k]] = k;
    }

    const Mt: number[][] = []; // scratch memory
    // Apply permutation. We make M into P*M*P^T
    for (let a = 0; a < M.length; ++a) {
      let ar = Pinv[M[a][0]];
      let ac = Pinv[M[a][1]];

      // we only store the upper-diagonal elements(since we assume matrix is symmetric, we only need to store these)
      // if permuted element is below diagonal, we simply transpose it.
      if (ac < ar) {
        const t = ac;
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
  const Ap = new Array<number>(n + 1);
  const Ai = new Array<number>(M.length);
  const Ax = new Array<number>(M.length);

  // count number of non-zero elements in columns.
  const LNZ: number[] = [];
  for (let i = 0; i < n; ++i) {
    LNZ[i] = 0;
  }
  for (const entry of M) {
    LNZ[entry[1]]++;
  }

  Ap[0] = 0;
  for (let i = 0; i < n; ++i) {
    Ap[i + 1] = Ap[i] + LNZ[i];
  }

  const coloffset: number[] = [];
  for (let a = 0; a < n; ++a) {
    coloffset[a] = 0;
  }

  // go through all elements in M, and add them to sparse matrix A.
  for (const e of M) {
    const col = e[1];

    const adr = Ap[col] + coloffset[col];
    Ai[adr] = e[0];
    Ax[adr] = e[2];

    coloffset[col]++;
  }

  const D = new Array<number>(n);
  const Y = new Array<number>(n);
  const Lp = new Array<number>(n + 1);
  const Parent = new Array<number>(n);
  const Lnz = new Array<number>(n);
  const Flag = new Array<number>(n);
  const Pattern = new Array<number>(n);
  const bp1 = new Array<number>(n);
  const x = new Array<number>(n);

  ldlSymbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag);

  const Lx = new Array<number>(Lp[n]);
  const Li = new Array<number>(Lp[n]);

  const d = ldlNumeric(
    n,
    Ap,
    Ai,
    Ax,
    Lp,
    Parent,
    Lnz,
    Li,
    Lx,
    D,
    Y,
    Pattern,
    Flag,
  );

  if (d === n) {
    return (b: NumberArray): number[] => {
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

export default prepare;
