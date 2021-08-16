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

@inline function exp_ps_m7_rho6_4(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_ps_m7_rho6_4!(A_copy)
end

@inline function exp_ps_m7_rho6_4!(A)
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
    # Computation order: B2 B3 B4 B5 Ba8_6 Bb6 B6 Bb7 B7 Bb8 B8 T2k10
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing B3 with operation: mult
    mul!(memslots3,memslots2,memslots1)
    # Computing B4 with operation: mult
    mul!(memslots4,memslots3,memslots1)
    # Computing B5 with operation: mult
    mul!(memslots5,memslots4,memslots1)
    # Computing Ba8_6 = x*I+x*B5
    coeff1=0.40308943770549005
    coeff2=1.0
    memslots6 .= coeff2.*memslots5
    mul!(memslots6,true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=7.647163731819816e-13
    coeff2=4.779477332387385e-14
    coeff3=2.8114572543455206e-15
    coeff4=1.5619206968586225e-16
    coeff5=8.22063524662433e-18
    coeff6=4.110317623312165e-19
    memslots7 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4 .+ coeff6.*memslots5
    mul!(memslots7,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots8,memslots5,memslots7)
    # Deallocating Bb6 in slot 7
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B6
    coeff1=2.755731922398589e-7
    coeff2=2.505210838544172e-8
    coeff3=2.08767569878681e-9
    coeff4=1.6059043836821613e-10
    coeff5=1.1470745597729725e-11
    coeff6=1.0
    # Smart lincomb recycle B6
    memslots8 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4 .+ coeff6.*memslots8
    mul!(memslots8,true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots7,memslots5,memslots8)
    # Deallocating B5 in slot 5
    # Deallocating Bb7 in slot 8
    # Computing Bb8 = x*I+x*A+x*B2+x*B3+x*B4+x*B7
    coeff1=0.008333333333333333
    coeff2=0.001388888888888889
    coeff3=0.0001984126984126984
    coeff4=2.48015873015873e-5
    coeff5=2.7557319223985893e-6
    coeff6=1.0
    # Smart lincomb recycle B7
    memslots7 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4 .+ coeff6.*memslots7
    mul!(memslots7,true,I*coeff1,true,true)
    # Computing B8 with operation: mult
    mul!(memslots5,memslots6,memslots7)
    # Deallocating Ba8_6 in slot 6
    # Deallocating Bb8 in slot 7
    # Computing T2k10 = x*I+x*A+x*B2+x*B3+x*B4+x*B8
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=0.16666666666666666
    coeff5=0.041666666666666664
    coeff6=1.0
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4 .+ coeff6.*memslots5
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 4
    # Deallocating B8 in slot 5
    return memslots1 # Returning T2k10
end

