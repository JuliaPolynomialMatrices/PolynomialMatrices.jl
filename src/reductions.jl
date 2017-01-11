# Computes the degree of each column of a polynomial matrix
function col_degree{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = degree(p)
  num_col = size(p,2)

  k = fill(-1,1,num_col)
  for i = max_deg:-1:0
    v = coeffs(p)[i]
    for j = 1:num_col
      if k[j] < 0 && !all(v[:,j] .== zero(T))
        k[j] = i
      end
    end
  end
  k
end

# Computes the degree of each row of a polynomial matrix
function row_degree{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = degree(p)
  num_row = size(p,1)

  k = fill(-1,num_row,1)
  for i = max_deg:-1:0
    v = coeffs(p)[i]
    for j = 1:num_row
      if k[j] < 0 && !all(v[j,:] .== zero(T))
        k[j] = i
      end
    end
  end
  k
end

# Computes the degree of each column, and the highest-column-degree
# coefficient matrix of a polynomial matrix
function high_col_deg_matrix{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = degree(p)
  num_col = size(p,2)

  k   = fill(-1,1,num_col)
  Phc = zeros(T,p.dims)
  for i = max_deg:-1:0
    c = coeffs(p)[i]
    for j = 1:num_col
      if k[j] < 0 && !all(c[:,j] .== zero(T))
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
  max_deg = degree(p)
  num_row = size(p,1)

  k   = fill(-1,num_row,1)
  Phr = zeros(p.dims)
  for i = max_deg:-1:0
    c = coeffs(p)[i]
    for j = 1:num_row
      if k[j] < 0 && !all(c[j,:] .== zero(T))
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
function colred{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  N < 2 || size(p,1) ≥ size(p,2) ||
    error("colred: Polynomial matrix is not full column rank")

  p_temp = copy(p)
  c       = p_temp.coeffs          # Dictionary of coefficient matrices of p
  num_col = N < 2 ? 1 : size(p,2)  # Number of columns of p
  U       = PolyMatrix(eye(T,num_col),p.var)

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
  size(p1,2) == size(p2,2) || (N1 < 2 && N2 < 2) ||
    error("colred: Both polynomial matrices should have the same number of columns")
  p1.var == p2.var ||
    error("colred: Both polynomial matrices should be in the same variable")
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
  (N < 2 && size(p,1) ≤ 1) || size(p,1) ≤ size(p,2) ||
    error("rowred: Polynomial matrix is not full row rank")
  p_temp  = copy(p)
  c       = coeffs(p_temp)  # Dictionary of coefficient matrices of p
  num_row = size(p,1)      # Number of rows of p
  U       = PolyMatrix(eye(T,num_row),p.var)

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
  size(p1,1) == size(p2,1) ||
    error("rowred: Both polynomial matrices should have the same number of rows")
  p1.var == p2.var ||
    error("rowred: Both polynomial matrices should be in the same variable")
  (N1 < 2 && size(p1,1) ≤ 1) || size(p1,1) ≤ size(p1,2) ||
    error("rowred: Polynomial matrix `p1` is not full row rank")

  p1_temp  = copy(p1)
  p2_temp  = copy(p2)
  c1       = coeffs(p_temp1)  # Dictionary of coefficient matrices of p1
  c2       = coeffs(p_temp2)  # Dictionary of coefficient matrices of p2
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
        c2[max_temp-k1[row]+i][Nmax,:] += n[row] / n[Nmax] * c2[i][row,:]
      end
    end

    # Reset collection indN
    fill!(indN, 0)
  end
end
