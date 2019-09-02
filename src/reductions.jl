"""
    gcrd(p₁, p₂, [iterative::Bool=true, dᵤ::Int]) -> R, V₁, V₂

Finds a greatest common right divisor of `p₁` `p₂`.
The method returns unimodular `V₁` and `V₂` such that `p₁ = V₁ R`, and
`p₂ = V₂ R`, where `R` is upper triangular.

By default the iterative [version 3'][1] is used to produce the
triangularization. If `iterative` is to `false`, the [method 3][1] tries to
triangularize `p` with a reduction matrix of degree `dᵤ`. By default `dᵤ` is set
high enough to guarantee triangularization if it is possible.

# Examples
```julia
julia> s = variable("s");
julia> p₁ = PolyMatrix([s^2 zero(s); zero(s) s^2]);
julia> p₂ = PolyMatrix([one(s) s+1]);
julia> R, V₁, V₂ = gcrd(p₁,p₂);
julia> R
2x2 PolyArray{Float64,2}:
  Poly(-1.0)  Poly(-1.0 - 1.0⋅s)
  Poly(-0.0)  Poly(-1.1547⋅s^2)
```

# References

-  [1]: D. Henrion, M. Sebek "Reliable Numerical Methods for Polynomial Matrix
        Triangularization" IEEE Transactions on Automatic Control, vol. 44,
        no. 3, Mar. 1999.
"""
function gcrd(p₁::PolyMatrix{T1,M1,Val{W},N},
  p₂::PolyMatrix{T2,M2,Val{W},N}, iterative::Bool=true, dᵤ::Int=-1) where {T1,M1,W,N,T2,M2}
  n₁,m₁ = size(p₁)
  n₂,m₂ = size(p₂)
  m₁ == m₂ || throw(DimensionMismatch("the two polynomial matrices do not have the same number of columns"))
  R,U   = rtriang([p₁; p₂], iterative, dᵤ)
  detU, adjU = inv(U)
  V     = adjU/detU(0)
  V₁    = V[1:n₁,1:m₁]
  V₂    = V[n₁+(1:n₂),1:m₂]
  return R[1:m₁,1:m₁], V₁, V₂
end

"""
    gcld(p₁, p₂, [iterative::Bool=true, dᵤ::Int]) -> L, V₁, V₂

Finds a greatest common left divisor of `p₁` `p₂`.
The method returns unimodular `V₁` and `V₂` such that `p₁ = L V₁`, and
`p₂ = L V₂`, where `L` is lower triangular.

By default the iterative [version 3'][1] is used to produce the
triangularization. If `iterative` is to `false`, the [method 3][1] tries to
triangularize `p` with a reduction matrix of degree `dᵤ`. By default `dᵤ` is set
high enough to guarantee triangularization if it is possible.

# Examples
```julia
julia> s = variable("s");
julia> p₁ = PolyMatrix([s^2 zero(s); zero(s) s^2]);
julia> p₂ = PolyMatrix([one(s); s+1]);
julia> L, V₁, V₂ = gcld(p₁,p₂);
julia> L
2x2 PolyArray{Float64,2}:
  Poly(-1.0)          Poly(-0.0)
  Poly(-1.0 - 1.0⋅s)  Poly(-1.1547⋅s^2)
```

# References

-  [1]: D. Henrion, M. Sebek "Reliable Numerical Methods for Polynomial Matrix
        Triangularization" IEEE Transactions on Automatic Control, vol. 44,
        no. 3, Mar. 1999.
"""
function gcld(p₁::PolyMatrix{T1,M1,Val{W},N},
  p₂::PolyMatrix{T2,M2,Val{W},N}, iterative::Bool=true, dᵤ::Int=-1) where {T1,M1,W,N,T2,M2}
  n₁,m₁ = size(p₁)
  n₂,m₂ = size(p₂)
  n₁ == n₂ || throw(DimensionMismatch("the two polynomial matrices do not have the same number of columns"))
  L,U   = ltriang([p₁ p₂], iterative, dᵤ)
  detU, adjU = inv(U)
  V     = adjU/detU(0)
  V₁    = V[1:n₁,1:m₁]
  V₂    = V[1:n₁,m₁+(1:m₂)]
  return L[1:n₁,1:n₁], V₁, V₂
end

"""
    hermite(p, iterative::Bool=true, dᵤ) -> H, U

Hermite form of a Polynomial Matrix.
The method returns unimodular `H` and `U` such that `p U = H`, where `H` is
in Hermite form of `p`.

By default `p` is triangularized using the iterative [version 3'][1] (see `ltriang`).
The triangular form is subsequently brought into Hermite form using the Euclidian
algorithm.

# Examples
```julia
julia> s = variable("s");
julia> p = PolyMatrix([-s^3-2s^2+1 -(s+1)^2; (s+2)^2*(s+1) zero(s)]);
julia> H,U = hermite(p);
julia> H
2x2 PolyArray{Float64,2}:
  Poly(1.0 + 1.0⋅s)                      Poly(0.0)
  Poly(4.0 + 8.0⋅s + 5.0⋅s^2 + 1.0⋅s^3)  Poly(4.0 + 12.0⋅s + 13.0⋅s^2 + 6.0⋅s^3 + 4.0⋅s^4)
```

# References

-  [1]: D. Henrion, M. Sebek "Reliable Numerical Methods for Polynomial Matrix
        Triangularization" IEEE Transactions on Automatic Control, vol. 44,
        no. 3, Mar. 1999.
"""
function hermite(p::PolyMatrix{T1,M,Val{W},N}, iterative::Bool=true, dᵤ::Int=-1) where {T1,M,W,N}
  L,U,d = _ltriang(p, iterative, dᵤ)
  n,m   = size(p)

  # scale diagonal elements first
  Σ = [findfirst(x -> x!=0, L[:,k]) for k in 1:m]
  U1 = Diagonal([1 ./ L[Σ[k],k] for k in 1:m])
  U = U*U1
  L = L*U1

  # reduce order
  U2 = Matrix{Float64}(I,m,m)
  for i in 2:m
    for j in 1:i-1
      σ = Σ[i]
      U2[i,j] = -L[σ,j]
    end
  end

  U = U*U2
  L = L*U2

  Lᵣ = _unshift(L,d)
  return PolyMatrix(Lᵣ, (n,m), Val{W}; reverse=true), PolyMatrix(U, (m,m), Val{W}; reverse=true)
end

"""
    ltriang(p, iterative::Bool=true, dᵤ) -> L, U

Polynomial Matrix triangularization based on Sylvester matrix [method 3][1].
The method returns unimodular `L` and `U` such that `p U = L`, where `L` is lower triangular.

By default the iterative [version 3'][1] is used. If `iterative` is to `false`,
the [method 3][1] tries to triangularize `p` with a reduction matrix of degree `dᵤ`.
By default `dᵤ` is set high enough to guarantee triangularization if it is possible.

Note that not necessarily the hermite form is returned (see `hermite`).

# Examples
```julia
julia> s = variable("s");
julia> p = PolyMatrix([s-1 s^2-1; 2 2s+2; 0 3]);
julia> L,U = ltriang(p);
julia> L
3x1 PolyArray{Float64,2}:
  Poly(1.22474 - 1.22474⋅s)  Poly(0)
  Poly(-2.44949)             Poly(0)
  Poly(-1.22474)             Poly(1.73205)
```

# References

-  [1]: D. Henrion, M. Sebek "Reliable Numerical Methods for Polynomial Matrix
        Triangularization" IEEE Transactions on Automatic Control, vol. 44,
        no. 3, Mar. 1999.
"""
function ltriang(p::PolyMatrix{T1,M,Val{W},N}, iterative::Bool=true, dᵤ::Int=-1) where {T1,M,W,N}
  n,m = size(p)
  if n < m || rank(p) < m
    pₑ = vcat(p, PolyMatrix(Matrix{Float64}(undef,m,m), (m,m), Val{W}))
  else
    pₑ = p
  end
  L,U,d = _ltriang(pₑ, iterative, dᵤ)
  L = _unshift(L[1:n*(d+1),1:m],d)
  return PolyMatrix(L, (n,m), Val{W}; reverse=true), PolyMatrix(U, (m,m), Val{W}; reverse=true)
end

"""
    rtriang(p, iterative::Bool=true, dᵤ::Int) -> R, U

Polynomial Matrix triangularization based on Sylvester matrix [method 3][1].
The method returns unimodular `R` and `U` such that `U p = R`, where `R` is upper triangular.

By default the iterative [version 3'][1] is used. If `iterative` is to `false`,
the [method 3][1] tries to triangularize `p` with a reduction matrix of degree `dᵤ`.
By default `dᵤ` is set high enough to guarantee triangularization if it is possible.

Note that not necessarily the hermite form is returned (see `hermite`).

# Examples
```julia
julia> s = variable("s");
julia> p = PolyMatrix([s-1 2 0; s^2-1 2s+2 3]);
julia> R,U = rtriang(p);
julia> R
3x1 PolyArray{Float64,2}:
  Poly(1.22474 - 1.22474⋅s)  Poly(-2.44949)  Poly(-1.22474)
  Poly(0)                    Poly(0)         Poly(1.73205)
```

# References

-  [1]: D. Henrion, M. Sebek "Reliable Numerical Methods for Polynomial Matrix
        Triangularization" IEEE Transactions on Automatic Control, vol. 44,
        no. 3, Mar. 1999.
"""
function rtriang(p::PolyMatrix{T1,M,Val{W},N}, iterative::Bool=true, dᵤ::Int=-1) where {T1,M,W,N}
  L,U = ltriang(transpose(p), iterative, dᵤ)
  return transpose(L), transpose(U)
end

function _unshift(L::AbstractMatrix,d::Int)
  n,r = divrem(size(L,1), d+1)
  r == 0 || throw(DimensionMismatch())
  SL = zero(L)
  for i in 0:d
    for j in 0:n-1
      SL[(d-i)*n+j+1,:] = L[(j)*(d+1)+d-i+1,:]
    end
  end
  return SL
end

function _ltriang(p::PolyMatrix{T1,M,Val{W},N}, iterative::Bool=true, dᵤ::Int=-1) where {T1,M,W,N}
  # allow user defined dᵤ
  if !iterative && dᵤ < 0
    dᵤ = _mindegree(p)
  elseif dᵤ < 0
    dᵤ = 0
  end
  d   = degree(p)+dᵤ
  n,m = size(p)

  # tolerance for extracting ϵshape
  ϵ = sqrt(n)*sqrt(m)*d*1e-16
  T = float(T1)
  # construct row permuted sylvester matrix
  mₛ = m*(dᵤ+1)
  Rd = zeros(T, n*(d+1), mₛ)
  for (k,v) in coeffs(p)
    for i in 0:n-1, j in 0:dᵤ
      Rd[n*(d+1)-i*(d+1)-k-dᵤ+j, j*m .+ (1:m)] = v[n-i,:]
    end
  end
  # should be changed when support for 0.4 drops (lq not in 0.4)
  q,L = qr(transpose(Rd))
  L   = L'
  U   = q'

  triangularshape = false
  j = 1
  Σb = zeros(Int,mₛ)
  while !triangularshape && j < 2mₛ
    for i in 1:mₛ
      # TODO if index is zero p is not of full rank.
      # Should we call the method with Identity appended?
      Σb[i] = findfirst(x->abs(x) > ϵ, L[:,i])
      L[1:Σb[i]-1,i] = zeros(T, Σb[i]-1)
    end

    # sort in decreasing order
    v = sortperm(Σb)
    U2 = zeros(T, mₛ, mₛ)
    for k in 1:mₛ
      U2[k,v[k]] = one(T)
    end
    L = L*transpose(U2)
    U = U2*U
    Σb = Σb[v]

    # check triangular shape
    i = 2
    triangularshape = true
    while i <= mₛ
      if Σb[i] == Σb[i-1]
        qi,Li = qr(L[Σb[i]:end,i-1:end]')
        U[i-1:end,:]         = transpose(qi)*Up1  = Poly([1])
p2  = Poly([2,1,3])
p3  = Poly([2,3,4])
m   = [p1 p2; p2 p1]
pm1 = PolyMatrix(m, Val{:x})
#@inferred PolyMatrix(m, Val{:x})

d = Dict(0=>[1 2;2 1], 1=>[0 1;1 0], 2=>[0 3;3 0])
@test pm1 == PolyMatrix(d)
#@inferred PolyMatrix(d, Val{:x})
@test_throws DimensionMismatch PolyMatrix(Dict{Int,Matrix{Int}}())

d = Dict(0=>[1 2;2 1], 1=>[0 1;1 0; 0 0])
@test_throws DimensionMismatch PolyMatrix(d)

degreepm2 = 8
ny  = 2
nu  = 2
A   = randn(ny*(degreepm2+1),nu)
B   = Matrix{Float64}(I,2,2)
[i-1:end,:]
        L[Σb[i]:end,i-1:end] = Li'
        triangularshape = false
      end
      i += 1
    end
    j += 1
  end

  C = [zeros(Int,0) for i in 1:n]
  # calculate index sets
  for j in eachindex(Σb)
    σ    = Σb[j]
    it,r = divrem(σ, d+1)
    i    = r > 0 ? it+1 : it
    push!(C[i], j)
  end

  # The following line is suggested by testing and not by the paper ?!
  # It handles cases where the first rows in the
  # triangular form are zero
  Σb = Σb .- Σb[1] .+ 1

  Σ = zeros(Int,0)
  for i in 1:m*(dᵤ+1)
    if Σb[i] > n
      continue
    end
    if !isempty(C[Σb[i]])
      push!(Σ,i)
    end
  end

  Uₜ = transpose(U)
  if length(Σ) < m
    if iterative && dᵤ < _mindegree(p)
      return _ltriang(p, iterative, dᵤ+1)
    else
      throw(ErrorException("ltriang: failed to triangularize"))
    end
  else
    k = [maximum(C[i]) for i in Σ]
    Uᵣ = Uₜ[:,k]
    Lᵣ = L[:,k]
    # truncate accorting to ϵ
    for i in eachindex(Uᵣ)
      Uᵣ[i] = abs(Uᵣ[i]) > ϵ ? Uᵣ[i] : zero(T)
    end
    for i in eachindex(Lᵣ)
      Lᵣ[i] = abs(Lᵣ[i]) > ϵ ? Lᵣ[i] : zero(T)
    end
    return Lᵣ, Uᵣ, d
  end
end

function _mindegree(p::PolyMatrix)
  m = minimum(size(p))
  rowdegs = sort(vec(row_degree(p)))
  coldegs = sort(vec(col_degree(p)))
  return min(sum(rowdegs[end-m+2:end]), sum(coldegs[end-m+2:end]))
end

# Computes the degree of each column of a polynomial matrix
function col_degree(p::PolyMatrix{T,M,O,N}) where {T,M,O,N}
  max_deg = degree(p)
  num_col = size(p,2)
  c       = coeffs(p)

  k = fill(-1,1,num_col)
  for i = max_deg:-1:0
    if haskey(c,i)
      v = c[i]
      for j = 1:num_col
        if k[j] < 0 && !all(v[:,j] .== zero(T))
          k[j] = i
        end
      end
    end
  end
  k
end

# Computes the degree of each row of a polynomial matrix
function row_degree(p::PolyMatrix{T,M,O,N}) where {T,M,O,N}
  max_deg = degree(p)
  num_row = size(p,1)
  c       = coeffs(p)

  k = fill(-1,num_row,1)
  for i = max_deg:-1:0
    if haskey(c,i)
      v = c[i]
      for j = 1:num_row
        if k[j] < 0 && !all(v[j,:] .== zero(T))
          k[j] = i
        end
      end
    end
  end
  k
end

# Computes the degree of each column, and the highest-column-degree
# coefficient matrix of a polynomial matrix
function high_col_deg_matrix(p::PolyMatrix{T,M,O,N}) where {T,M,O,N}
  max_deg = degree(p)
  num_col = size(p,2)
  c       = coeffs(p)

  k   = fill(-1,1,num_col)
  Phc = zeros(T,size(p))
  for i = max_deg:-1:0
    if haskey(c,i)
      v = c[i]
      for j = 1:num_col
        if k[j] < 0 && !all(v[:,j] .== zero(T))
          k[j] = i
          Phc[:,j] = v[:,j]
        end
      end
    end
  end
  return k, Phc
end

# Computes the degree of each row, and the highest-row-degree
# coefficient matrix of a polynomial matrix
function high_row_deg_matrix(p::PolyMatrix{T,M,O,N}) where {T,M,O,N}
  max_deg = degree(p)
  num_row = size(p,1)
  c       = coeffs(p)

  k   = fill(-1,num_row,1)
  Phr = zeros(T,size(p))
  for i = max_deg:-1:0
    if haskey(c,i)
      v = c[i]
      for j = 1:num_row
        if k[j] < 0 && !all(v[j,:] .== zero(T))
          k[j] = i
          Phr[j,:] = v[j,:]
        end
      end
    end
  end
  return k, Phr
end

# Determines if a polynomial matrix is column proper (or "column reduced")
is_col_proper(p::PolyMatrix{T,M,O,1}) where {T,M,O} = true
function is_col_proper(p::PolyMatrix{T,M,O,N}) where {T,M,O,N}
  Phc = high_col_deg_matrix(p)[2]
  return rank(Phc) == size(p,2)
end

# Determines if a polynomial matrix is row proper (or "row reduced")
function is_row_proper(p::PolyMatrix)
  Phr = high_row_deg_matrix(p)[2]
  return rank(Phr) == size(p,1)
end

# Computes the column reduced form of a polynomial matrix
# (via Wolovich's method)
# NOTE: Wolovich's method is known to be numerically unstable (Geurts and Praagman, 1996);
# It would be preferable to implement Geurts-Praagman's method
# NOTE: Should the procedure end with an error if p is not full rank, or simply
# indicate this as an output argument (a la SLICOT)?
function colred(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}
  N < 2 || size(p,1) ≥ size(p,2) ||
    error("colred: Polynomial matrix is not full column rank")

  p_temp = copy(p)
  c       = coeffs(p_temp)         # Dictionary of coefficient matrices of p
  num_col = N < 2 ? 1 : size(p,2)  # Number of columns of p
  U       = PolyMatrix(Matrix{Float64}(undef,num_col,num_col), Val{W})

  indN    = zeros(Int,num_col)  # Collection of non-zero entries of n
  while true
    k, Phc = high_col_deg_matrix(p_temp)
    nPhc   = nullspace(Phc)
    if size(nPhc,2) == zero(T)
      return p_temp, U
    end
    n = view(nPhc,:,1)  # One vector from the nullspace of Phc

    # Find all nonzero entries of n, and among those,
    # the one associated with the column of highest degree
    num_nz = 0               # Number of elements in indN
    Nmax   = 0               # Index of the column of p with highest degree
    max_temp = -1            # Maximum degree found so far
    for i = 1:num_col
      if n[i] != 0
        num_nz += 1
        indN[num_nz] = i
        if k[i] > max_temp
          Nmax = i
          max_temp = k[i]
        end
      end
    end

    # If there are less than 2 nonzero entries in indN,
    # the polynomial matrix is not full rank (because it has a column of zeros)
    if num_nz < 2
      error("colred: Polynomial matrix is not full column rank")
    end

    # Unimodular matrix Utemp
    Utemp = SortedDict(0 => Matrix{Float64}(num_col,num_col))
    for i = 1:max_temp-minimum(k[indN[1:num_nz]])
      insert!(Utemp, i, zeros(T,num_col,num_col))
    end

    # Perform column reduction
    for j = 1:num_nz
      if j == Nmax
        continue
      end
      col = indN[j]

      # Update coefficient matrices
      for i = 0:k[col]
        c[max_temp-k[col]+i][:,Nmax] += n[col] / n[Nmax] * c[i][:,col]
      end

      # Update Utemp
      Utemp[max_temp-k[col]][col,Nmax] = n[col] / n[Nmax]
    end

    # Update unimodular transformation matrix U
    U = U*PolyMatrix(Utemp, (num_col,num_col), Val{W})

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes the simultaneous column reduced form of two polynomial matrices `p1` and `p2`,
# according to `p1`(via Wolovich's method)
function colred(p1::PolyMatrix{T,M1,Val{W},N1},
  p2::PolyMatrix{T,M2,Val{W},N2}) where {T,M1,M2,W,N1,N2}
  size(p1,2) == size(p2,2) || (N1 < 2 && N2 < 2) ||
    error("colred: Both polynomial matrices should have the same number of columns")
  N1 < 2 || size(p1,1) ≥ size(p1,2) ||
    error("colred: Polynomial matrix `p1` is not full column rank")

  p1_temp  = copy(p1)
  p2_temp  = copy(p2)
  c1       = coeffs(p1_temp)          # Dictionary of coefficient matrices of p1
  c2       = coeffs(p2_temp)          # Dictionary of coefficient matrices of p2
  num_col  = N1 < 2 ? 1 : size(p1,2)  # Number of columns of p1 and p2

  indN    = zeros(Int,num_col)  # Collection of non-zero entries of n
  while true
    k1, Phc = high_col_deg_matrix(p1_temp)
    k2      = col_degree(p2_temp)
    nPhc    = nullspace(Phc)
    if size(nPhc,2) == zero(T)
      return p1_temp, p2_temp
    end

    n = view(nPhc,:,1)  # One vector from the nullspace of Phc

    # Find all nonzero entries of n, and among those,
    # the one associated with the column of highest degree
    num_nz   = 0             # Number of elements in indN
    Nmax     = 0             # Index of the column of p1 with highest degree
    max_temp = -1            # Maximum degree found so far
    for i = 1:num_col
      if n[i] != 0
        num_nz += 1
        indN[num_nz] = i
        if k1[i] > max_temp
          Nmax     = i
          max_temp = k1[i]
        end
      end
    end

    # If there are less than 2 nonzero entries in indN,
    # the polynomial matrix is not full rank (because it has a column of zeros)
    if num_nz < 2
      error("colred: Polynomial matrix is not full rank")
    end

    # Perform column reduction
    for j = 1:num_nz
      if j == Nmax
        continue
      end
      col = indN[j]

      # Update coefficient matrices of p1
      for i = 0:k1[col]
        c1[max_temp-k1[col]+i][:,Nmax] += n[col] / n[Nmax] * c1[i][:,col]
      end

      # Update coefficient matrices of p2
      for i = 0:k2[col]
        if haskey(c2,i)
          new_key = max_temp-k1[col]+i
          new_update = n[col] / n[Nmax] * c2[i][:,col]
          if haskey(c2,new_key)
            c2[new_key][:,Nmax] += new_update
          else
            v2 = zeros(first(c2)[2])
            v2[:,Nmax] += new_update
            insert!(c2, new_key, v2)
          end
        end
      end
    end

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes the row reduced form of a polynomial matrix
# (via Wolovich's method)
function rowred(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}
  (N < 2 && size(p,1) ≤ 1) || size(p,1) ≤ size(p,2) ||
    error("rowred: Polynomial matrix is not full row rank")
  p_temp  = copy(p)
  c       = coeffs(p_temp)  # Dictionary of coefficient matrices of p
  num_row = size(p,1)      # Number of rows of p
  U       = PolyMatrix(Matrix{Float64}(num_row,num_row), Val{W})

  indN    = zeros(Int,num_row)  # Collection of non-zero entries of n
  while true
    k, Phr = high_row_deg_matrix(p_temp)
    nPhr   = nullspace(Phr')
    if size(nPhr,2) == zero(T)
      return p_temp, U
    end

    n = view(nPhr,:,1)  # One vector from the nullspace of Phc

    # Find all nonzero entries of n, and among those,
    # the one associated with the row of highest degree
    num_nz = 0               # Number of elements in indN
    Nmax   = 0               # Index of the row of p with highest degree
    max_temp = -1            # Maximum degree found so far
    for i = 1:num_row
      if n[i] != 0
        num_nz += 1
        indN[num_nz] = i
        if k[i] > max_temp
          Nmax = i
          max_temp = k[i]
        end
      end
    end

    # If there are less than 2 nonzero entries in indN,
    # the polynomial matrix is not full rank (because it has a row of zeros)
    if num_nz < 2
      error("rowred: Polynomial matrix is not full rank")
    end

    # Unimodular matrix Utemp
    Utemp = SortedDict(0 => Matrix{Float64}(num_row,num_row))
    #insert!(Utemp, 0, )
    for i = 1:max_temp-minimum(k[indN[1:num_nz]])
      insert!(Utemp, i, zeros(T,num_row,num_row))
    end

    # Perform row reduction
    for j = 1:num_nz
      if j == Nmax
        continue
      end
      row = indN[j]

      # Update coefficient matrices
      for i = 0:k[row]
        c[max_temp-k[row]+i][Nmax,:] += n[row] / n[Nmax] * c[i][row,:]
      end

      # Update Utemp
      Utemp[max_temp-k[row]][Nmax,row] = n[row] / n[Nmax]
    end

    # Update unimodular transformation matrix U
    U = PolyMatrix(Utemp, (num_row,num_row), Val{W})*U

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes the simultaneous row reduced form of two polynomial matrices `p1` and `p2`,
# according to `p1`(via Wolovich's method)
function rowred(p1::PolyMatrix{T,M1,Val{W},N1},
  p2::PolyMatrix{T,M2,Val{W},N2}) where {T,M1,M2,W,N1,N2}
  size(p1,1) == size(p2,1) ||
    error("rowred: Both polynomial matrices should have the same number of rows")
  (N1 < 2 && size(p1,1) ≤ 1) || size(p1,1) ≤ size(p1,2) ||
    error("rowred: Polynomial matrix `p1` is not full row rank")

  p1_temp  = copy(p1)
  p2_temp  = copy(p2)
  c1       = coeffs(p1_temp)  # Dictionary of coefficient matrices of p1
  c2       = coeffs(p2_temp)  # Dictionary of coefficient matrices of p2
  num_row  = size(p1,1)      # Number of rows of p1 and p2

  indN    = zeros(Int,num_row)  # Collection of non-zero entries of n
  while true
    k1, Phr = high_row_deg_matrix(p1_temp)
    k2      = row_degree(p2_temp)
    nPhr    = nullspace(Phr')
    if size(nPhr,2) == zero(T)
      return p1_temp, p2_temp
    end

    n = view(nPhr,:,1)  # One vector from the nullspace of Phr

    # Find all nonzero entries of n, and among those,
    # the one associated with the row of highest degree
    num_nz   = 0             # Number of elements in indN
    Nmax     = 0             # Index of the row of p1 with highest degree
    max_temp = -1            # Maximum degree found so far
    for i = 1:num_row
      if n[i] != 0
        num_nz += 1
        indN[num_nz] = i
        if k1[i] > max_temp
          Nmax     = i
          max_temp = k1[i]
        end
      end
    end

    # If there are less than 2 nonzero entries in indN,
    # the polynomial matrix is not full rank (because it has a row of zeros)
    if num_nz < 2
      error("rowred: Polynomial matrix is not full rank")
    end

    # Perform row reduction
    for j = 1:num_nz
      if j == Nmax
        continue
      end
      row = indN[j]

      # Update coefficient matrices of p1
      for i = 0:k1[row]
        c1[max_temp-k1[row]+i][Nmax,:] += n[row] / n[Nmax] * c1[i][row,:]
      end

      # Update coefficient matrices of p2
      for i = 0:k2[row]
        if haskey(c2,i)
          new_key = max_temp-k1[row]+i
          new_update = n[row] / n[Nmax] * c2[i][row,:]
          if haskey(c2,new_key)
            c2[new_key][Nmax,:] += new_update
          else
            v2 = zeros(first(c2)[2])
            v2[Nmax,:] += new_update
            insert!(c2, new_key, v2)
          end
        end
      end
    end

    # Reset collection indN
    fill!(indN, 0)
  end
end
