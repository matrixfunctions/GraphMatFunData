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

@inline function exp_ps_rho2_7(A)
    return exp_ps_rho2_7!(copy(A))
end

@inline function exp_ps_rho2_7!(A)
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
    # Computation order: B2 B3 B4 Bb5 B5 Bb6 B6 Bb7 Ba7_5 B7 T2k9
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing B3 with operation: mult
    mul!(memslots[3],memslots[2],memslots[1])
    # Computing B4 with operation: mult
    mul!(memslots[4],memslots[3],memslots[1])
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=2.08767569878681e-9
    coeff2=1.6059043836821613e-10
    coeff3=1.1470745597729725e-11
    coeff4=7.647163731819816e-13
    coeff5=4.779477332387385e-14
    memslots[5] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4]
    mul!(memslots[5],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[6],memslots[4],memslots[5])
    # Deallocating Bb5 in slot 5
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B5
    coeff1=2.48015873015873e-5
    coeff2=2.7557319223985893e-6
    coeff3=2.755731922398589e-7
    coeff4=2.505210838544172e-8
    coeff5=1.0
    # Smart lincomb recycle B5
    memslots[6] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[6]
    mul!(memslots[6],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[5],memslots[4],memslots[6])
    # Deallocating Bb6 in slot 6
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B6
    coeff1=0.041666666666666664
    coeff2=0.008333333333333333
    coeff3=0.001388888888888889
    coeff4=0.0001984126984126984
    coeff5=1.0
    # Smart lincomb recycle B6
    memslots[5] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[5]
    mul!(memslots[5],true,I*coeff1,true,true)
    # Computing Ba7_5 = x*I+x*B4
    coeff1=0.20154471885274503
    coeff2=1.0
    # Smart lincomb recycle B4
    memslots[4] .= coeff2.*memslots[4]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots[6],memslots[4],memslots[5])
    # Deallocating Ba7_5 in slot 4
    # Deallocating Bb7 in slot 5
    # Computing T2k9 = x*I+x*A+x*B2+x*B3+x*B7
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=0.16666666666666666
    coeff5=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[6]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B7 in slot 6
    return memslots[1] # Returning T2k9
end

