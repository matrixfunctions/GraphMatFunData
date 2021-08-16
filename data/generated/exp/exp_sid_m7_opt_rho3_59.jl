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

@inline function exp_sid_m7_opt_rho3_59(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_sid_m7_opt_rho3_59!(A_copy)
end

@inline function exp_sid_m7_opt_rho3_59!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=10
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
    memslots10 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 Ba7 Bb7 B7 Ba8 Bb8 B8 T2k10
    # Computing Bb2 = x*I+x*A
    coeff1=-0.0017737608141714808
    coeff2=1.005235674329508
    memslots2 .= coeff2.*memslots1
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=-0.0017737608141714808
    coeff2=1.005235674329508
    memslots3 .= coeff2.*memslots1
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=0.0005326111625306977
    coeff2=1.0012299933973516
    coeff3=0.006925033919227824
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=-0.00016401351885017156
    coeff2=0.0005433816370422478
    coeff3=1.0012270195002069
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots5,memslots2,memslots3)
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=9.544368792763562e-5
    coeff2=0.9999194856569542
    coeff3=0.0036194298213211165
    coeff4=-0.00033449839479958844
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=-1.1883273082175028e-5
    coeff2=3.080294633483493e-5
    coeff3=9.623360783744299e-5
    coeff4=0.9998936828701804
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots6,memslots2,memslots3)
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-6.0508143407070435e-6
    coeff2=1.0000102174714505
    coeff3=0.0028500766374399016
    coeff4=-0.00013226806180853073
    coeff5=1.0208196110243714e-5
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-4.786848660402673e-7
    coeff2=8.582321314936397e-8
    coeff3=5.9167378426502e-6
    coeff4=-5.991737523333378e-6
    coeff5=0.999999999087183
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots7,memslots2,memslots3)
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-1.5738212245123645e-8
    coeff2=3.740920280371873e-8
    coeff3=5.267964654623867e-7
    coeff4=-3.46635995191657e-6
    coeff5=3.708627877122971e-6
    coeff6=0.9999999999895876
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-6.288899284895006e-6
    coeff2=1.1355716740441329e-6
    coeff3=5.567958735815622e-8
    coeff4=2.487015589947048e-9
    coeff5=1.0836272945889375e-10
    coeff6=1.1368496406363211e-11
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots8,memslots2,memslots3)
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.00022214082957652806
    coeff2=5.947941137388186
    coeff3=0.4220821708139979
    coeff4=0.024845741420535115
    coeff5=0.0011101347889166462
    coeff6=3.168781522592016e-5
    coeff7=1.003633457505924
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-0.00014732216151790324
    coeff2=-0.004469043351396654
    coeff3=0.0326061354999038
    coeff4=0.0028514008243872304
    coeff5=0.0001864934813936358
    coeff6=3.852437481436439e-5
    coeff7=0.9846670594297366
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots9,memslots2,memslots3)
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 3
    # Computing Ba8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.2015451143023047
    coeff2=11.106844258571114
    coeff3=2.9912660184614275
    coeff4=0.2891134043244349
    coeff5=0.02649378029184827
    coeff6=0.0027421911866186664
    coeff7=61.7597688836338
    coeff8=0.9846397447218478
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.00028017073414700555
    coeff2=-0.0013783167325222296
    coeff3=0.007816923517808423
    coeff4=0.0011284750735624438
    coeff5=8.405278388325319e-5
    coeff6=1.3621609668499039e-5
    coeff7=0.9977252838919913
    coeff8=3.4108017105321994e-5
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B8 with operation: mult
    mul!(memslots10,memslots2,memslots3)
    # Deallocating Ba8 in slot 2
    # Deallocating Bb8 in slot 3
    # Computing T2k10 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=0.9999420628667928
    coeff2=0.9990039337658048
    coeff3=0.5079542977446111
    coeff4=0.08249596102472202
    coeff5=0.005411384044110833
    coeff6=7.400879751137793e-5
    coeff7=-4.7284436864949066e-5
    coeff8=-0.0023105536908983168
    coeff9=0.9977227168431884
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9 .+ coeff9.*memslots10
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    # Deallocating B7 in slot 9
    # Deallocating B8 in slot 10
    return memslots1 # Returning T2k10
end

