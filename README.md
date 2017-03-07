# PolynomialMatrices

[![Unix][unix-img]][unix-link]
[![Windows][win-img]][win-link]
[![Coveralls][ca-img]][ca-link]
[![Codecov][cc-img]][cc-link]

[unix-img]: https://img.shields.io/travis/neveritt/PolynomialMatrices.jl/master.svg?label=unix
[unix-link]: https://travis-ci.org/neveritt/PolynomialMatrices.jl
[win-img]: https://img.shields.io/appveyor/ci/neveritt/PolynomialMatrices-jl/master.svg?label=windows
[win-link]: https://ci.appveyor.com/project/neveritt/PolynomialMatrices-jl/branch/master
[ca-img]: https://img.shields.io/coveralls/neveritt/PolynomialMatrices.jl/master.svg?label=coveralls
[ca-link]: https://coveralls.io/github/neveritt/PolynomialMatrices.jl?branch=master
[cc-img]: https://img.shields.io/codecov/c/github/neveritt/PolynomialMatrices.jl/master.svg?label=codecov
[cc-link]: https://codecov.io/gh/neveritt/PolynomialMatrices.jl?branch=master

This package is a work in progress to provide univariate polynomial matrix calculations.

### Description

Some of the implemented and yet to be implemented functionality

- [x] Core PolynomialMatrix
  - [x] basic PolynomialMatrix interface
  - [x] use SortedDict for internal representation
- [ ] Math operations
  - [x] +
  - [x] -
  - [x] *
  - [ ] /
- [x] Indexing:
  - [x] getindex
  - [x] setindex!
  - [x] Iteration
- [x] Constructor
  - [x] constructor from matrix of Polynomials.Poly
  - [x] constructor form dict of AbstractMatrix
  - [x] constructor from tall AbstractMatrix
  - [x] constructor from 3-dimensional AbstractMatrix
- [x] Functions
  - [x] degree, rank, determinant, transpose, inv
  - [x] row and column reduction
  - [ ] gcd
- [ ] transformations
  - [x] triangularization
  - [x] transformation to hermite form
  - [ ] transformation to smith form
