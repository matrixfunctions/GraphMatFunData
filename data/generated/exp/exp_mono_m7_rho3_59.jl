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

@inline function exp_mono_m7_rho3_59(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_mono_m7_rho3_59!(A_copy)
end

@inline function exp_mono_m7_rho3_59!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=8
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    memslots7 = similar(A,T)
    memslots8 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: B2 B3 B4 B5 B6 B7 B8 T2k10
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing B3 with operation: mult
    mul!(memslots3,memslots2,memslots1)
    # Computing B4 with operation: mult
    mul!(memslots4,memslots3,memslots1)
    # Computing B5 with operation: mult
    mul!(memslots5,memslots4,memslots1)
    # Computing B6 with operation: mult
    mul!(memslots6,memslots5,memslots1)
    # Computing B7 with operation: mult
    mul!(memslots7,memslots6,memslots1)
    # Computing B8 with operation: mult
    mul!(memslots8,memslots7,memslots1)
    # Computing T2k10 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=0.16666666666666666
    coeff5=0.041666666666666664
    coeff6=0.008333333333333333
    coeff7=0.001388888888888889
    coeff8=0.0001984126984126984
    coeff9=2.48015873015873e-5
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4 .+ coeff6.*memslots5 .+ coeff7.*memslots6 .+ coeff8.*memslots7 .+ coeff9.*memslots8
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 4
    # Deallocating B5 in slot 5
    # Deallocating B6 in slot 6
    # Deallocating B7 in slot 7
    # Deallocating B8 in slot 8
    return memslots1 # Returning T2k10
end

