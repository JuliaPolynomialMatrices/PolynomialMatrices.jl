# Computes the degree of each column of a polynomial matrix
function col_degree{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = order(p)
  num_col = N < 2 ? 1 : p.dims[2]

  k   = fill(-1,1,num_col)
  for i = max_deg:-1:0
    c = p.coeffs[i]
    for j = 1:num_col
      if k[j] < 0 && !all(c[:,j] .== 0)
        k[j] = i
      end
    end
  end
  return k
end

# Computes the degree of each row of a polynomial matrix
function row_degree{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = order(p)
  num_row = p.dims[1]

  k   = fill(-1,num_row,1)
  for i = max_deg:-1:0
    c = p.coeffs[i]
    for j = 1:num_row
      if k[j] < 0 && !all(c[j,:] .== 0)
        k[j] = i
      end
    end
  end
  return k
end

# Computes the degree of each column, and the highest-column-degree
# coefficient matrix of a polynomial matrix
function high_col_deg_matrix{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = order(p)
  num_col = N < 2 ? 1 : p.dims[2]

  k   = fill(-1,1,num_col)
  Phc = zeros(p.dims)
  for i = max_deg:-1:0
    c = p.coeffs[i]
    for j = 1:num_col
      if k[j] < 0 && !all(c[:,j] .== 0)
        k[j] = i
        Phc[:,j] = c[:,j]
      end
    end
  end
  return k, Phc
end

# Computes the degree of each row, and the highest-row-degree
# coefficient matrix of a polynomial matrix
function high_row_deg_matrix{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = order(p)
  num_row = p.dims[1]

  k   = fill(-1,num_row,1)
  Phr = zeros(p.dims)
  for i = max_deg:-1:0
    c = p.coeffs[i]
    for j = 1:num_row
      if k[j] < 0 && !all(c[j,:] .== 0)
        k[j] = i
        Phr[j,:] = c[j,:]
      end
    end
  end
  return k, Phr
end

# Determines if a polynomial matrix is column proper (or "column reduced")
is_col_proper{T,M,O}(p::PolyMatrix{T,M,O,1}) = true
function is_col_proper{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  k, Phc = high_col_deg_matrix(p)
  return rank(Phc) == p.dims[2]
end

# Determines if a polynomial matrix is row proper (or "row reduced")
function is_row_proper(p::PolyMatrix)
  k, Phr = high_row_deg_matrix(p)
  return rank(Phr) == p.dims[1]
end

# Computes the column reduced form of a polynomial matrix
# (via Wolovich's method)
# NOTE: Wolovich's method is known to be numerically unstable (Geurts and Praagman, 1996);
# It would be preferable to implement Geurts-Praagman's method
# NOTE: Should the procedure end with an error if p is not full rank, or simply
# indicate this as an output argument (a la SLICOT)?
function colred{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  @assert N < 2 || p.dims[1] ≥ p.dims[2]
    "colred: Polynomial matrix is not full column rank"

  p_temp = copy(p)
  c       = p_temp.coeffs          # Dictionary of coefficient matrices of p
  num_col = N < 2 ? 1 : p.dims[2]  # Number of columns of p
  U       = PolyMatrix(eye(T,num_col),p.var)

  indN    = zeros(Int,num_col)  # Collection of non-zero entries of n
  while true
    k, Phc = high_col_deg_matrix(p_temp)
    nPhc   = nullspace(Phc)
    if size(nPhc,2) == 0
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
    Utemp = SortedDict(Dict{Int,AbstractMatrix{}}())
    insert!(Utemp, 0, eye(T,num_col))
    for i = 1:max_temp-minimum(k[indN])
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
    U = U*PolyMatrix(Utemp,(num_col,num_col),p.var)

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes the simultaneous column reduced form of two polynomial matrices `p1` and `p2`,
# according to `p1`(via Wolovich's method)
function colred{T,M1,M2,O1,O2,N1,N2}(p1::PolyMatrix{T,M1,O1,N1},
  p2::PolyMatrix{T,M2,O2,N2})
  @assert p1.dims[2] == p2.dims[2] || (N1 < 2 && N2 < 2)
    "colred: Both polynomial matrices should have the same number of columns"
  @assert p1.var == p2.var
    "colred: Both polynomial matrices should be in the same variable"
  @assert N1 < 2 || p1.dims[1] ≥ p1.dims[2]
    "colred: Polynomial matrix `p1` is not full column rank"

  p1_temp  = copy(p1)
  p2_temp  = copy(p2)
  c1       = p1_temp.coeffs           # Dictionary of coefficient matrices of p1
  c2       = p2_temp.coeffs           # Dictionary of coefficient matrices of p2
  num_col  = N1 < 2 ? 1 : p1.dims[2]  # Number of columns of p1 and p2

  indN    = zeros(Int,num_col)  # Collection of non-zero entries of n
  while true
    k1, Phc = high_col_deg_matrix(p1_temp)
    k2      = col_degree(p2_temp)
    nPhc    = nullspace(Phc)
    if size(nPhc,2) == 0
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
        c2[max_temp-k1[col]+i][:,Nmax] += n[col] / n[Nmax] * c2[i][:,col]
      end
    end

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes the row reduced form of a polynomial matrix
# (via Wolovich's method)
function rowred{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  @assert (N < 2 && p.dims[1] ≦ 1) || p.dims[1] ≦ p.dims[2]
    "rowred: Polynomial matrix is not full row rank"
  p_temp  = copy(p)
  c       = p_temp.coeffs  # Dictionary of coefficient matrices of p
  num_row = p.dims[1]      # Number of rows of p
  U       = PolyMatrix(eye(num_row),p.var)

  indN    = zeros(Int,num_row)  # Collection of non-zero entries of n
  while true
    k, Phr = high_row_deg_matrix(p_temp)
    nPhr   = nullspace(Phr')
    if size(nPhr,2) == 0
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
    Utemp = SortedDict(Dict{Int,AbstractMatrix{}}())
    insert!(Utemp, 0, eye(T,num_row))
    for i = 1:max_temp-minimum(k[indN])
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
    U = PolyMatrix(Utemp,(num_row,num_row),p.var)*U

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes the simultaneous row reduced form of two polynomial matrices `p1` and `p2`,
# according to `p1`(via Wolovich's method)
function rowred{T,M1,M2,O1,O2,N1,N2}(p1::PolyMatrix{T,M1,O1,N1},
  p2::PolyMatrix{T,M2,O2,N2})
  @assert p1.dims[1] == p2.dims[1]
    "rowred: Both polynomial matrices should have the same number of rows"
  @assert p1.var == p2.var
    "rowred: Both polynomial matrices should be in the same variable"
  @assert (N1 < 2 && p1.dims[1] ≤ 1) || p1.dims[1] ≤ p1.dims[2]
    "rowred: Polynomial matrix `p1` is not full row rank"

  p1_temp  = copy(p1)
  p2_temp  = copy(p2)
  c1       = p1_temp.coeffs  # Dictionary of coefficient matrices of p1
  c2       = p2_temp.coeffs  # Dictionary of coefficient matrices of p2
  num_row  = p1.dims[1]      # Number of rows of p1 and p2

  indN    = zeros(Int,num_row)  # Collection of non-zero entries of n
  while true
    k1, Phr = high_row_deg_matrix(p1_temp)
    k2      = row_degree(p2_temp)
    nPhr    = nullspace(Phr')
    if size(nPhr,2) == 0
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
        c2[max_temp-k1[row]+i][Nmax,:] += n[row] / n[Nmax] * c2[i][row,:]
      end
    end

    # Reset collection indN
    fill!(indN, 0)
  end
end
