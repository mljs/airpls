'use strict';

const airPLS = require('.src/index');

var a = new Array(10).fill(1);

for (let i = 2; i < 3; i++) {
    a[i] += 5;
}

let result = airPLS(a, {maxIterations: 4, lambda: 100});
console.log(result.baseline);
console.log(result.error)