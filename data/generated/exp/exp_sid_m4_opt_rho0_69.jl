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

@inline function exp_sid_m4_opt_rho0_69(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_sid_m4_opt_rho0_69!(A_copy)
end

@inline function exp_sid_m4_opt_rho0_69!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=7
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    memslots7 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 T2k7
    # Computing Bb2 = x*I+x*A
    coeff1=9.94823321550182e-5
    coeff2=1.0000210350390846
    memslots2 .= coeff2.*memslots1
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=9.94823321550182e-5
    coeff2=1.0000210350390846
    memslots3 .= coeff2.*memslots1
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=8.745327870027237e-6
    coeff2=-4.197092420395855e-8
    coeff3=0.9999999958392435
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=-1.0660766661003063e-5
    coeff2=0.0029440745098486642
    coeff3=0.0004021729193984898
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots5,memslots2,memslots3)
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=4.022059101053733e-5
    coeff2=0.4019191862379914
    coeff3=-0.00873832228961145
    coeff4=1.0000192435512933
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=2.1672609863472577e-5
    coeff2=6.626362532960593e-5
    coeff3=0.03233860555846545
    coeff4=1.0000837252012718
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots6,memslots2,memslots3)
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.00011854240934339541
    coeff2=2.2246677805188155
    coeff3=0.2615074307076249
    coeff4=5.768958229431865
    coeff5=0.9990123769097652
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-7.86695052350628e-5
    coeff2=-0.04148569168401734
    coeff3=0.023305174055381638
    coeff4=-0.023059406187239697
    coeff5=0.9998635489305508
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots7,memslots2,memslots3)
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing T2k7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.9999999940884431
    coeff2=0.9999715234150552
    coeff3=0.5918961310449276
    coeff4=-3.2734161536753197
    coeff5=10.408133461599245
    coeff6=0.9998619329256802
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    return memslots1 # Returning T2k7
end

