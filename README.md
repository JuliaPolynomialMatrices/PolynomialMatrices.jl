# PolynomialMatrices

This package is a work in progress to provide univariate polynomial matrix calculations.

### Build Status and Code Coverage

# -  Build status: [![Build Status][build-img]][build-link]
# -  Code coverage: [![Coveralls][ca-img]][ca-link] [![Codecov][cc-img]][cc-link]

[build-img]:  https://travis-ci.org/neveritt/PolynomialMatrices.jl.svg?branch=master
[build-link]: https://travis-ci.org/neveritt/PolynomialMatrices.jl
[ca-img]: https://coveralls.io/repos/github/neveritt/PolynomialMatrices.jl/badge.svg?branch=master
[ca-link]: https://coveralls.io/github/neveritt/PolynomialMatrices.jl?branch=master
[cc-img]: https://codecov.io/gh/neveritt/PolynomialMatrices.jl/branch/master/graph/badge.svg
[cc-link]: https://codecov.io/gh/neveritt/PolynomialMatrices.jl

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

