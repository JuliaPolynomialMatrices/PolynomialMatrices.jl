# PolynomialMatrices

Short description

### Build Status and Code Coverage

# -  Build status: [![Build Status][build-img]][build-link]
# -  Code coverage: [![Coveralls][ca-img]][ca-link] [![Codecov][cc-img]][cc-link]

[build-img]:  https://travis-ci.org/KTH-AC/PolynomialMatrices.jl.svg?branch=master
[build-link]: https://travis-ci.org/KTH-AC/PolynomialMatrices.jl
[ca-img]: https://coveralls.io/repos/github/KTH-AC/PolynomialMatrices.jl/badge.svg?branch=master
[ca-link]: https://coveralls.io/github/KTH-AC/PolynomialMatrices.jl?branch=master
[cc-img]: https://codecov.io/gh/KTH-AC/PolynomialMatrices.jl/branch/master/graph/badge.svg
[cc-link]: https://codecov.io/gh/KTH-AC/PolynomialMatrices.jl

### Description

#### TODO's

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
- [x] Functions
	- [x] transpose
  - [x] ctranspose
	- [x] matrix functions (inv, transpose, etc...) (could be more)
