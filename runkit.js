'use strict';

const airPLS = require('ml-airpls');

let y = [1, 1, 1, 1, 3, 6, 3, 1, 1, 1];
let x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

let result = airPLS(x, y, {maxIterations: 10, lambda: 100});

console.log(result.baseline);
