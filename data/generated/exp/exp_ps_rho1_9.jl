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

@inline function exp_ps_rho1_9(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_ps_rho1_9!(A_copy)
end

@inline function exp_ps_rho1_9!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=6
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: B2 B3 B4 Bb5 B5 Bb6 B6 T2k8
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing B3 with operation: mult
    mul!(memslots3,memslots2,memslots1)
    # Computing B4 with operation: mult
    mul!(memslots4,memslots3,memslots1)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=2.48015873015873e-5
    coeff2=2.7557319223985893e-6
    coeff3=2.755731922398589e-7
    coeff4=2.505210838544172e-8
    coeff5=2.08767569878681e-9
    memslots5 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4
    mul!(memslots5,true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots6,memslots4,memslots5)
    # Deallocating Bb5 in slot 5
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B5
    coeff1=0.041666666666666664
    coeff2=0.008333333333333333
    coeff3=0.001388888888888889
    coeff4=0.0001984126984126984
    coeff5=1.0
    # Smart lincomb recycle B5
    memslots6 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots6
    mul!(memslots6,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots5,memslots4,memslots6)
    # Deallocating B4 in slot 4
    # Deallocating Bb6 in slot 6
    # Computing T2k8 = x*I+x*A+x*B2+x*B3+x*B6
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=0.16666666666666666
    coeff5=1.0
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots5
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B6 in slot 5
    return memslots1 # Returning T2k8
end

