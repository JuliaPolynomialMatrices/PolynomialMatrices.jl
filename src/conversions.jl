promote_rule{T1,T2,M1,M2,V,N}(::Type{PolyMatrix{T1,M1,Val{V},N}},
  ::Type{PolyMatrix{T2,M2,Val{V},N}}) =
  PolyMatrix{promote_type(T1, T2), promote_type(M1,M2), Val{V}, N}

function _convert{T1,N,T2,M1,M2,V}(::Type{PolyMatrix{T1,M1,Val{V},N}},
  p::PolyMatrix{T2,M2,Val{V},N})
  r = PolyMatrix( SortedDict(Dict{Int,M1}()), size(p), Val{V})
  for (k,c) in coeffs(p)
    r.coeffs[k] = map(x->convert(T1,x),c)
  end
  r
end

@generated function convert{T1,N,T2,M1,M2,V}(::Type{PolyMatrix{T1,M1,Val{V},N}},
  p::PolyMatrix{T2,M2,Val{V},N})
  if !(M1 <: AbstractArray{T1,N})
    return :(error("convert: first two parameters of first argument are incompatible"))
  elseif T1 == T2 && M1 == M2
    return :(p)
  else
    return :(_convert(PolyMatrix{T1,M1,Val{V},N}, p::PolyMatrix{T2,M2,Val{V},N}))
  end
end
