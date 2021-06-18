using LinearAlgebra

struct ValueOne; end
ValueOne()

# Compute X <- a X + b I.
function matfun_axpby!(X,a,b,Y::UniformScaling)
    m,n=size(X)
    if ~(a isa ValueOne)
        rmul!(X,a)
    end
    @inbounds for i=1:n
        X[i,i]+=(b isa ValueOne) ? 1 : b
    end
end

@inline function exp_ps_m2(A)
    return exp_ps_m2!(copy(A))
end

@inline function exp_ps_m2!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=4
    memslots=Vector{Matrix{T}}(undef,max_memslots)
    n=size(A,1)
    for j=1:max_memslots
        memslots[j]=Matrix{T}(undef,n,n)
    end
    # The first slots are precomputed nodes [:A]
    memslots[1]=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: A2 P1 C0 P0
    # Computing A2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing P1 = x*I+x*A+x*A2
    coeff1=0.5
    coeff2=0.16666666666666666
    coeff3=0.041666666666666664
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[2]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing C0 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating P1 in slot 3
    # Deallocating A2 in slot 2
    # Computing P0 = x*I+x*A+x*C0
    coeff1=1.0
    coeff2=1.0
    coeff3=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating C0 in slot 4
    return memslots[1] # Returning P0
end

