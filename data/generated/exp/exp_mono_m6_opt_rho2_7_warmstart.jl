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

@inline function exp_mono_m6_opt_rho2_7_warmstart(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_mono_m6_opt_rho2_7_warmstart!(A_copy)
end

@inline function exp_mono_m6_opt_rho2_7_warmstart!(A)
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
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 Ba7 Bb7 B7 T2k9
    # Computing Bb2 = x*I+x*A
    coeff1=-0.00059913427279684
    coeff2=1.0025622563474483
    memslots2 .= coeff2.*memslots1
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=-0.00059913427279684
    coeff2=1.0025622563474483
    memslots3 .= coeff2.*memslots1
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=-0.0001186386176465649
    coeff2=-0.0007786622491555297
    coeff3=0.9997089384902429
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=-0.0008389425622496766
    coeff2=1.001103922131545
    coeff3=0.016675233270506518
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots5,memslots2,memslots3)
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-1.0901226987553663e-5
    coeff2=-0.00024610324529555493
    coeff3=-0.0012675065756547724
    coeff4=0.9998216665236417
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.001913396449826322
    coeff2=1.0003404217204235
    coeff3=0.03779228316890699
    coeff4=0.0010026121023330907
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots6,memslots2,memslots3)
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-1.0845198429404047e-6
    coeff2=-3.927444922940479e-5
    coeff3=-0.0001726204582222337
    coeff4=-0.0018914370462158503
    coeff5=0.9999810493229588
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.002428107252826335
    coeff2=0.9991866421557342
    coeff3=0.06048518858862958
    coeff4=0.0030848564056926066
    coeff5=6.530788756304574e-5
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots7,memslots2,memslots3)
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=5.845729699785629e-7
    coeff2=-4.804333361952689e-6
    coeff3=-2.174097745805705e-5
    coeff4=-0.00017747138030956434
    coeff5=-0.0023957202877644673
    coeff6=0.9999986701463371
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.0011640774645410574
    coeff2=0.9979163876428353
    coeff3=0.07557668582102052
    coeff4=0.00513929879938074
    coeff5=0.00016385434755408757
    coeff6=6.796434505014407e-8
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots8,memslots2,memslots3)
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=1.2510851844021513e-7
    coeff2=-1.0549528366134809e-7
    coeff3=-2.488531283091839e-6
    coeff4=-4.148575000794752e-5
    coeff5=4.326460110080258e-5
    coeff6=0.0011668444495142097
    coeff7=0.9999999802111392
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-8.09922010751223e-8
    coeff2=0.9970626583923993
    coeff3=0.09545225884934809
    coeff4=0.008678577298423432
    coeff5=0.0004953140491424132
    coeff6=1.3254090277403565e-5
    coeff7=3.942451324846202e-8
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots9,memslots2,memslots3)
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 3
    # Computing T2k9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.999999803981589
    coeff2=1.0006179487544273
    coeff3=0.4979281278966493
    coeff4=0.16592558923333828
    coeff5=0.03871151836039848
    coeff6=0.0061833894649526416
    coeff7=0.0006066317390421412
    coeff8=3.0487323320460782e-5
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    # Deallocating B7 in slot 9
    return memslots1 # Returning T2k9
end

