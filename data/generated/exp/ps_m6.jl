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

@inline function ps_m6(A)
    return ps_m6!(copy(A))
end

@inline function ps_m6!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=6
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
    # Computation order: A2 A3 A4 P3 C2 P2 C1 P1 C0 P0
    # Computing A2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing A3 with operation: mult
    mul!(memslots[3],memslots[1],memslots[2])
    # Computing A4 with operation: mult
    mul!(memslots[4],memslots[1],memslots[3])
    # Computing P3 = x*I+x*A+x*A2+x*A3+x*A4
    coeff1=2.08767569878681e-9
    coeff2=1.6059043836821613e-10
    coeff3=1.1470745597729725e-11
    coeff4=7.647163731819816e-13
    coeff5=4.779477332387385e-14
    memslots[5] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4]
    mul!(memslots[5],true,I*coeff1,true,true)
    # Computing C2 with operation: mult
    mul!(memslots[6],memslots[5],memslots[4])
    # Deallocating P3 in slot 5
    # Computing P2 = x*I+x*A+x*A2+x*A3+x*C2
    coeff1=2.48015873015873e-5
    coeff2=2.7557319223985893e-6
    coeff3=2.755731922398589e-7
    coeff4=2.505210838544172e-8
    coeff5=1.0
    # Smart lincomb recycle C2
    memslots[6] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[6]
    mul!(memslots[6],true,I*coeff1,true,true)
    # Computing C1 with operation: mult
    mul!(memslots[5],memslots[6],memslots[4])
    # Deallocating P2 in slot 6
    # Computing P1 = x*I+x*A+x*A2+x*A3+x*C1
    coeff1=0.041666666666666664
    coeff2=0.008333333333333333
    coeff3=0.001388888888888889
    coeff4=0.0001984126984126984
    coeff5=1.0
    # Smart lincomb recycle C1
    memslots[5] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[5]
    mul!(memslots[5],true,I*coeff1,true,true)
    # Computing C0 with operation: mult
    mul!(memslots[6],memslots[5],memslots[4])
    # Deallocating P1 in slot 5
    # Deallocating A4 in slot 4
    # Computing P0 = x*I+x*A+x*A2+x*A3+x*C0
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=0.16666666666666666
    coeff5=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[6]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating A2 in slot 2
    # Deallocating A3 in slot 3
    # Deallocating C0 in slot 6
    return memslots[1] # Returning P0
end

