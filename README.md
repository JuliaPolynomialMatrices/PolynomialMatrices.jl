# PolynomialMatrices

[![Build Status](https://github.com/JuliaPolynomialMatrices/PolynomialMatrices.jl/workflows/CI/badge.svg)](https://github.com/JuliaPolynomialMatrices/PolynomialMatrices.jl/actions)
[![Codecov](https://codecov.io/gh/JuliaPolynomialMatrices/PolynomialMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaPolynomialMatrices/PolynomialMatrices.jl)

The `PolynomialMatrices` Julia package provides a support for calculations with univariate [polynomial matrices](https://en.wikipedia.org/wiki/Polynomial_matrix), that is, matrices whose entries are univariate polynomials, such as

```math
        [ 3s²+2s+1     1 ]
P(s) =  |                |
        [    2s        s ].
```

Tighly related to systems of linear ordinary differential or difference equations, these mathematical objects are useful in disciplines such as automatic control and signal processing.

## Polynomial matrix in Julia as an `Array` of polynomials

A matrix of univariate polynomials can be created in Julia by combining the functionality of the standard Julia `Array` (or `Matrix`) type and the `Poly` type provided by the [`Polynomials`](https://github.com/JuliaMath/Polynomials.jl) package. Our example polynomial matrix can be entered in Julia as

```julia
julia> using Polynomials

julia> M = [Poly([1, 2, 3]) Poly([1]); Poly([0,2]) Poly([0,1])]
2×2 Array{Polynomials.Poly{Int64,2}:
  Poly(1 + 2⋅x + 3⋅x^2)  Poly(1)
  Poly(2⋅x)              Poly(x)
```

The trouble with the resulting matrix is that the number of operations we can do over such objest is quite limited. If, for example, we want to check column-reducedness or perform triangularization, plain Julia has no functionality for that. And now `PolynomialMatrices` package enters the stage...

## Polynomial matrix in Julia as `PolyMatrix`

Using `PolynomialMatrices` package, we convert our array of polynomials into a `PolyMatrix` type as in

```julia
julia> using PolynomialMatrices

julia> P = PolyMatrix(M)
2×2 PolyArray{Int64,2}:
  Poly(1 + 2⋅x + 3⋅x^2)  Poly(1)
  Poly(2⋅x)              Poly(x)
```

for which many operations are now defined. For example, the matrix can be evaluated at a given value of its independent variable

```julia
julia> P(1)
2×2 Array{Int64,2}:
  6  1
  2  1
```

Or, a polynomial matrix `P` can be transformed into an upper triangular `R` by premultplication by unimodular `U` using

```julia
julia> R,U = rtriang(P)
(Poly{Float64}[Poly(-0.5454545454545455) Poly(-0.5454545454545454 + 0.6666666666666669*x + 0.5757575757575755*x^2 - 0.36363636363636365*x^3); Poly(0.0) Poly(-0.2357022603955159*x + 0.4714045207910319*x^2 + 0.7071067811865472*x^3)], Poly{Float64}[Poly(-0.5454545454545454 + 0.24242424242424268*x) Poly(0.4242424242424243 + 0.5757575757575756*x - 0.36363636363636365*x^2); Poly(-0.4714045207910317*x) Poly(0.23570226039551567 + 0.47140452079103196*x + 0.7071067811865472*x^2)])

julia> R
2×2 PolyMatrix{Float64,Array{Float64,2},Val{:x},2}:
 Poly(-0.545455)  Poly(-0.545455 + 0.666667*x + 0.575758*x^2 - 0.363636*x^3)
 Poly(0.0)        Poly(-0.235702*x + 0.471405*x^2 + 0.707107*x^3)
```
And some more computation with polynomial matrices is offered by the package.

## Polynomial matrix viewed (and entered) as a matrix polynomial
A very useful interpretation of a polynomial matrix is that of a matrix polynomial. That is, a polynomial whose coefficients are not just numbers but matrices. Our original example can thus be written as

```math
        [ 1  1 ]   [ 2  0 ]     [ 3  0 ]
P(s) =  |      | + |      | s + |      | s²
        [ 0  0 ]   [ 2  1 ]     [ 0  0 ].
```

Hence, a natural way to enter a polynomial matrix in Julia is by entering a 3D array of coefficients matrices (of the corresponding matrix polynomial).

```julia
julia> A = zeros(Int,2,2,3);
julia> A[:,:,1] = [1 1;0 0];
julia> A[:,:,2] = [2 0;2 1];
julia> A[:,:,3] = [3 0;0 0];

julia> P = PolyMatrix(A,:s)
2×2 PolyArray{Int64,2}:
  Poly(1 + 2⋅s + 3⋅s^2)  Poly(1)
  Poly(2⋅s)              Poly(s)
```

## Polynomial matrix stored internally (and entered) as a dictionary

A `PolyMatrix` is implemented and stored as a `dict`, mapping from the powers (of the variables) to the coefficient matrices. It can thus also be constructed from a `dict` as in

```julia
julia> d = Dict(0=>[1 1;0 0], 1=>[2 0;2 1], 2=>[3 0;0 0]);
julia> P = PolyMatrix(d,:s)
2×2 PolyArray{Int64,2}:
  Poly(1 + 2⋅s + 3⋅s^2)  Poly(1)
  Poly(2⋅s)              Poly(s)
```
Individual coefficient matrices can be accessed accordingly --- the coefficient dictionary is obtained using `coeffs` function and the individual coefficient matrices are accessed using keys. For example, the coefficient matrix with the 1st power of the variable can be obtained using

```julia
julia> C = coeffs(P)
DataStructures.SortedDict{Int64,Array{Int64,2},Base.Order.ForwardOrdering} with 3 entries:
  0 => [1 1; 0 0]
  1 => [2 0; 2 1]
  2 => [3 0; 0 0]

julia> C[1]
2×2 Array{Int64,2}:
 2  0
 2  1
```

## List of functions for `PolyMatrix` objects

The functions for polynomial matrices implemented in `PolynomialMatrices` package are:

### Inquiry about parameters of the polynomial matrix
* `coeffs`: dictionary of coefficient matrices, keys are the powers
* `degree`, `col_degree`, `row_degree`: degree, column and row degrees
* `variable`, `vartype`: polynomial corresponding to the variable, symbol of the variable
* `high_col_deg_matrix`, `high_row_deg_matrix`: coefficient matrices corresponding to leading column and row degrees, respectively.

### Analysis
* `is_col_proper`, `is_row_proper`: checking if the matrix is column- and row-proper (also column- and row-reduced)

### Reductions, conversions
* `colred`, `rowred`: column and row degree reduction of a polynomial matrix
* `ltriang`, `rtriang`: conversion to a lower left and uppper right triangular polynomial matrix
* `hermite`: conversion to hermite form.

## Future plans
* `det` for computing the determinant of a polynomial matrix
* `roots` for computing the roots (or zeros) of a polynomial matrix
* some state space realization from a fraction of two polynomial matrices
* ...
* separate (but related) packages `PolynomialMatrixEquations` and `PolynomialMatrixFactorizations` are planned.

## `PolynomialMatrices` package is restricted to univariate polynomials only, for multivariate polynomials look elsewhere

As it is implemented now, `PolyMatrix` objects do not allow for mixing
different variables --- a `PolyMatrix` object can only operate together
with `PolyMatrix` objects with the same variable. For multivariate polynomials, you may want to check [`MultivariatePolynomials`](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) package.
