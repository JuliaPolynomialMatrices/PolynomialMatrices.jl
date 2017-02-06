promote_rule{T1,T2,M2,O,N}(::Type{PolyMatrix{T1,Array{T1,N},O,N}},
  ::Type{PolyMatrix{T2,M2,O,N}}) =
  PolyMatrix{promote_type(T1, T2), promote_type(M1,M2), O, N}

function _convert{T1,N,T2,M1,M2,O}(::Type{PolyMatrix{T1,M1,O,N}},
  p::PolyMatrix{T2,M2,O,N})
  r = PolyMatrix( SortedDict(Dict{Int,AbstractArray{T1,N}}()), size(p), p.var)
  for (k,c) in coeffs(p)
    r.coeffs[k] = map(x->convert(T1,x),c)
  end
  r
end

@generated function convert{T1,N,T2,M1,M2,O}(::Type{PolyMatrix{T1,M1,O,N}},
  p::PolyMatrix{T2,M2,O,N})
  if !(M1 <: AbstractArray{T1,N})
    return :(error("convert: first two parameters of first argument are incompatible"))
  elseif T1 == T2 && M1 == M2
    return :(p)
  else
    return :(_convert(PolyMatrix{T1,M1,O,N}, p::PolyMatrix{T2,M2,O,N}))
  end
end
