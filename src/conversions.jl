# PolyMatrix and PolyMatrix
promote_rule{T1,T2,M1,M2,W,N}(::Type{PolyMatrix{T1,M1,Val{W},N}},
  ::Type{PolyMatrix{T2,M2,Val{W},N}}) =
  PolyMatrix{promote_type(T1, T2), promote_type(M1,M2), Val{W}, N}

function _convert{T1,N,T2,M1,M2,W}(::Type{PolyMatrix{T1,M1,Val{W},N}},
  p::PolyMatrix{T2,M2,Val{W},N})
  cr = SortedDict(Dict{Int,M1}())
  for (k,c) in coeffs(p)
    insert!(cr, k, map(x->convert(T1,x),c))
  end
  PolyMatrix(cr, size(p), Val{W})
end

@generated function convert{T1,N,T2,M1,M2,W}(::Type{PolyMatrix{T1,M1,Val{W},N}},
  p::PolyMatrix{T2,M2,Val{W},N})
  if !(M1 <: AbstractArray{T1,N})
    return :(error("convert: first two parameters of first argument are incompatible"))
  elseif T1 == T2 && M1 == M2
    return :(p)
  else
    return :(_convert(PolyMatrix{T1,M1,Val{W},N}, p::PolyMatrix{T2,M2,Val{W},N}))
  end
end

# PolyMatrix and AbstractArray
promote_rule{T1,T2,M1,W,N}(::Type{PolyMatrix{T1,M1,Val{W},N}},
  ::Type{AbstractArray{T2,N}}) = PolyMatrix{promote_type(T1, T2), M1, Val{W}, N}

convert{T1,T2,M1,W,N}(::Type{PolyMatrix{T1,M1,Val{W},N}},
  p::AbstractArray{T2,N}) = PolyMatrx(p, size(p), Val{W})
