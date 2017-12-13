'use strict';

const airPLS = require('ml-airPLS');

var a = new Array(10).fill(1);

for (let i = 2; i < 4; i++) {
    a[i] += 5;
}

let result = airPLS(a, {maxIterations: 10, lambda: 100});

console.log(result.baseline);
