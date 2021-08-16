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


# Compute X <- a X + b Y.
function matfun_axpby!(X,a,b,Y)
    m,n=size(X)
    if ~(a isa ValueOne)
        rmul!(X,a)
    end
    @inbounds for i=1:m
        @inbounds for j=1:n
            if (b isa ValueOne)
                X[i,j]+=Y[i,j]
            else
                X[i,j]+=b*Y[i,j]
            end
        end
    end
end

@inline function exp_sastre_m8(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_sastre_m8!(A_copy)
end

@inline function exp_sastre_m8!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=9
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    memslots7 = similar(A,T)
    memslots8 = similar(A,T)
    memslots9 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: B2 B3 B4 B5 Bb6 B6 Ba7 Bb7 B7 Ba8 B8 Ba9 B9 T2k11
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing B3 with operation: mult
    mul!(memslots3,memslots1,memslots2)
    # Computing B4 with operation: mult
    mul!(memslots4,memslots1,memslots3)
    # Computing B5 with operation: mult
    mul!(memslots5,memslots1,memslots4)
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-1.023660713518307e-11
    coeff2=-4.508311519886735e-13
    coeff3=-1.980157255925737e-14
    coeff4=-9.210033748491798e-16
    coeff5=-6.140022498994532e-17
    memslots6 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4 .+ coeff5.*memslots5
    # Computing B6 with operation: mult
    mul!(memslots7,memslots5,memslots6)
    # Deallocating Bb6 in slot 6
    # Computing Ba7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-5.893435534477677e-5
    coeff2=-3.013961104055248e-6
    coeff3=-1.502070379373464e-7
    coeff4=-6.770221628797445e-9
    coeff5=-1.227011356117036e-10
    coeff6=1.0
    memslots6 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4 .+ coeff5.*memslots5 .+ coeff6.*memslots7
    # Computing Bb7 = x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-5.100472475630675e-7
    coeff2=-4.032817333361947e-8
    coeff3=-2.785084196756015e-9
    coeff4=-3.294026127901678e-10
    coeff5=1.0
    memslots8 .= coeff1.*memslots2 .+ coeff2.*memslots3 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots7
    # Computing B7 with operation: mult
    mul!(memslots9,memslots6,memslots8)
    # Deallocating Ba7 in slot 6
    # Deallocating Bb7 in slot 8
    # Computing Ba8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=2.755731922398589e-7
    coeff2=2.505210838544172e-8
    coeff3=2.08767569878681e-9
    coeff4=1.30531132637709e-10
    coeff5=7.55676813469492e-12
    coeff6=4.024189993755686e-13
    coeff7=-0.001023463999572971
    coeff8=1.0
    # Smart lincomb recycle B6
    memslots7 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4 .+ coeff6.*memslots5 .+ coeff7.*memslots7 .+ coeff8.*memslots9
    mul!(memslots7,true,I*coeff1,true,true)
    # Deallocating B7 in slot 9
    # Computing B8 with operation: mult
    mul!(memslots6,memslots7,memslots5)
    # Deallocating Ba8 in slot 7
    # Computing Ba9 = x*I+x*A+x*B2+x*B3+x*B4+x*B8
    coeff1=0.008333333333333333
    coeff2=0.001388888888888889
    coeff3=0.0001984126984126984
    coeff4=2.48015873015873e-5
    coeff5=2.7557319223985893e-6
    coeff6=1.0
    # Smart lincomb recycle B8
    memslots6 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4 .+ coeff6.*memslots6
    mul!(memslots6,true,I*coeff1,true,true)
    # Computing B9 with operation: mult
    mul!(memslots7,memslots6,memslots5)
    # Deallocating Ba9 in slot 6
    # Deallocating B5 in slot 5
    # Computing T2k11 = x*I+x*A+x*B2+x*B3+x*B4+x*B9
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=0.16666666666666666
    coeff5=0.041666666666666664
    coeff6=1.0
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4 .+ coeff6.*memslots7
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 4
    # Deallocating B9 in slot 7
    return memslots1 # Returning T2k11
end

