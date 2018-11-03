# airPLS

  [![NPM version][npm-image]][npm-url]
  [![build status][travis-image]][travis-url]
  [![Test coverage][codecov-image]][codecov-url]
  [![David deps][david-image]][david-url]
  [![npm download][download-image]][download-url]

Baseline correction using adaptive iteratively reweighted penalized least squares

It is an javascript implementation of [airPLS](https://github.com/zmzhang/airPLS/blob/master/airPLS_manuscript.pdf) using cholesky decomposition and reverse Cuthill-Mckee method for reducing the bandwidth of sparse linear systems, obtaining a fast baseline fitter. 

## Installation

`$ npm install ml-airPLS`

## [API Documentation](https://mljs.github.io/airPLS/)

## Example

```js
const airPls = require('ml-airPLS');

var yData = [0,0,0,2,5,2,0,0];
var {baseline, corrected, iteration, error} = airPLS(yData);
```

## License

[MIT](./LICENSE)

[npm-image]: https://img.shields.io/npm/v/ml-airPLS.svg?style=flat-square
[npm-url]: https://www.npmjs.com/package/ml-airPLS
[travis-image]: https://img.shields.io/travis/mljs/airPLS/master.svg?style=flat-square
[travis-url]: https://travis-ci.org/mljs/airPLS
[codecov-image]: https://img.shields.io/codecov/c/github/mljs/airPLS.svg?style=flat-square
[codecov-url]: https://codecov.io/gh/mljs/airPLS
[david-image]: https://img.shields.io/david/mljs/airPLS.svg?style=flat-square
[david-url]: https://david-dm.org/mljs/airPLS
[download-image]: https://img.shields.io/npm/dm/ml-airPLS.svg?style=flat-square
[download-url]: https://www.npmjs.com/package/ml-airPLS
