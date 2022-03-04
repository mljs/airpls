# airpls

  [![NPM version][npm-image]][npm-url]
  [![Test coverage][codecov-image]][codecov-url]
  [![npm download][download-image]][download-url]

Baseline correction using adaptive iteratively reweighted penalized least squares

It is an javascript implementation of [airpls](https://github.com/zmzhang/airPLS/blob/master/airPLS_manuscript.pdf) using cholesky decomposition and reverse Cuthill-Mckee method for reducing the bandwidth of sparse linear systems, obtaining a fast baseline fitter. 

## Installation

`$ npm install ml-airpls`

## [API Documentation](https://mljs.github.io/airpls/)

## Example

```js
const airpls = require('ml-airpls');

let y = [1, 1, 1, 1, 3, 6, 3, 1, 1, 1];
let x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

var {baseline, corrected, iteration, error} = airpls(x, y);
```

## License

[MIT](./LICENSE)

[npm-image]: https://img.shields.io/npm/v/ml-airpls.svg?style=flat-square
[npm-url]: https://www.npmjs.com/package/ml-airpls
[codecov-image]: https://img.shields.io/codecov/c/github/mljs/airpls.svg?style=flat-square
[codecov-url]: https://codecov.io/gh/mljs/airpls
[download-image]: https://img.shields.io/npm/dm/ml-airpls.svg?style=flat-square
[download-url]: https://www.npmjs.com/package/ml-airpls
