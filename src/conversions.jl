# PolyMatrix and PolyMatrix
promote_rule(::Type{PolyMatrix{T1,M1,Val{W},N}},
  ::Type{PolyMatrix{T2,M2,Val{W},N}}) where {T1,T2,M1,M2,W,N} =
  PolyMatrix{promote_type(T1, T2), promote_type(M1,M2), Val{W}, N}

function _convert(::Type{PolyMatrix{T1,M1,Val{W},N}},
  p::PolyMatrix{T2,M2,Val{W},N}) where {T1,N,T2,M1,M2,W}
  cr = SortedDict(Dict{Int,M1}())
  for (k,c) in coeffs(p)
    insert!(cr, k, map(x->convert(T1,x),c))
  end
  PolyMatrix(cr, size(p), Val{W})
end

@generated function convert(::Type{PolyMatrix{T1,M1,Val{W},N}},
  p::PolyMatrix{T2,M2,Val{W},N}) where {T1,N,T2,M1,M2,W}
  if !(M1 <: AbstractArray{T1,N})
    return :(error("convert: first two parameters of first argument are incompatible"))
  elseif T1 == T2 && M1 == M2
    return :(p)
  else
    return :(_convert(PolyMatrix{T1,M1,Val{W},N}, p::PolyMatrix{T2,M2,Val{W},N}))
  end
end

# PolyMatrix and AbstractArray
#promote_rule(::Type{PolyMatrix{T1,M1,Val{W},N}},
#  ::Type{AbstractArray{T2,N}}) where {T1,T2,M1,W,N} = PolyMatrix{promote_type(T1, T2), M1, Val{W}, N}

#convert(::Type{PolyMatrix{T1,M1,Val{W},N}},
#  p::AbstractArray{T2,N}) where {T1,T2,M1,W,N} = PolyMatrix(p, size(p), Val{W})
