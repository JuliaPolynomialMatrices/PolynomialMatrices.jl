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

## Description

`PolynomialMatrices` aims at supporting univariate polynomial matrix calculations.
The package extends `Polynomials`, and tries to provide a set of basic mathematical
functionality between `Number`s, `Poly`s and `PolyMatrix` objects.

## Basic Usage

The easiest way to construct a `PolyMatrix` object is to call its constructor
with a matrix of `Poly` objects. A `PolyMatrix` is implemented and stored as a
`dict`, mapping from the order to the coefficient matrices. It can thus also
be constructed from a `dict`.

As it is implemented now, `PolyMatrix` objects do not allow for mixing
different variables --- a `PolyMatrix` object's can only operate together
with `PolyMatrix` objects with the same variable.

For more information, check the documentation with `?PolynomialMatrices` command.

#### Example

```julia
julia> using Polynomials
julia> using PolynomialMatrices

julia> # construct PolyMatrix from matrix of polynomials
julia> m = [Poly([1, 2, 3]) Poly([1]); Poly([0,2]) Poly([0,1])]
2×2 Array{Polynomials.Poly{Int64,2}:
  Poly(1 + 2⋅x + 3⋅x^2)  Poly(1)
  Poly(2⋅x)              Poly(x)
julia> pm1 = PolyMatrix(m)
2×2 PolyArray{Int64,2}:
  Poly(1 + 2⋅x + 3⋅x^2)  Poly(1)
  Poly(2⋅x)              Poly(x)

julia> pm1(1)
2×2 Array{Int64,2}:
  6  1
  2  1

julia> # construct PolyMatrix from dictionary
julia> d = Dict(0=>[1 1;0 0], 1=>[2 0;2 1], 2=>[3 0;0 0]);
julia> pm2 = PolyMatrix(d,:s)
2×2 PolyArray{Int64,2}:
  Poly(1 + 2⋅s + 3⋅s^2)  Poly(1)
  Poly(2⋅s)              Poly(s)

julia> pm2(1)
2×2 Array{Int64,2}:
  6  1
  2  1

julia> # construct PolyMatrix from three-dimensional array
julia> a = zeros(Int,2,2,3);
julia> a[:,:,1] = [1 1;0 0];
julia> a[:,:,2] = [2 0;2 1];
julia> a[:,:,3] = [3 0;0 0];
julia> pm3 = PolyMatrix(a,:z)
2×2 PolyArray{Int64,2}:
  Poly(1 + 2⋅z + 3⋅z^2)  Poly(1)
  Poly(2⋅z)              Poly(z)

julia> pm3(1)
2×2 Array{Int64,2}:
  6  1
  2  1
```

### Convenience functions

Some functions exist for your convenience when working with `PolynomialMatrices`s.

Please read the corresponding documentation in Julia by issuing `?coeffs`, `?degree`,
`?variable`, `?vartype`, `?col_degree`, `?row_degree`, `?high_col_deg_matrix`,
`?high_row_deg_matrix`, `?is_col_proper`, `?is_row_proper`,
`?colred` or `?rowred`, `?ltriang`, `rtriang`, `hermite`.
