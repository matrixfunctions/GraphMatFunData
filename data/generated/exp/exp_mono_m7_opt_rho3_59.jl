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

@inline function exp_mono_m7_opt_rho3_59(A)
    return exp_mono_m7_opt_rho3_59!(copy(A))
end

@inline function exp_mono_m7_opt_rho3_59!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=10
    memslots=Vector{Matrix{T}}(undef,max_memslots)
    n=size(A,1)
    for j=1:max_memslots
        memslots[j]=Matrix{T}(undef,n,n)
    end
    # The first slots are precomputed nodes [:A]
    memslots[1]=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 Ba7 Bb7 B7 Ba8 Bb8 B8 T2k10
    # Computing Bb2 = x*I+x*A
    coeff1=1.3134712061518662e-5
    coeff2=1.002274766505785
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=1.3134712061518662e-5
    coeff2=1.002274766505785
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=-8.973017805624113e-5
    coeff2=0.00012005201731132445
    coeff3=0.9997367268186822
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=0.00012338976185732696
    coeff2=0.9996312949109102
    coeff3=0.018858289778215145
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots[5],memslots[2],memslots[3])
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-1.8416000490973626e-5
    coeff2=-8.104101755291844e-5
    coeff3=5.489518987562822e-5
    coeff4=0.9998500554706558
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=4.02314359582651e-6
    coeff2=0.9993459669804694
    coeff3=0.03609843676884375
    coeff4=0.0010584605114307261
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[2],memslots[3])
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-2.5675820144906785e-6
    coeff2=-1.618163016765615e-5
    coeff3=-3.4221784684677797e-6
    coeff4=2.4268115376905545e-5
    coeff5=0.9999837106332966
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=4.578521649969897e-5
    coeff2=1.0001351072704798
    coeff3=-0.013461797467855375
    coeff4=0.002324571658090968
    coeff5=0.0002315295921500602
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[2],memslots[3])
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-3.328521860384963e-7
    coeff2=-2.639951734327968e-6
    coeff3=-1.1160650661114392e-6
    coeff4=4.185994260757166e-6
    coeff5=7.151219728856269e-5
    coeff6=0.9999995656134203
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-0.00011293613093022973
    coeff2=0.9990923104443914
    coeff3=0.047438292970867935
    coeff4=0.0033538018355920127
    coeff5=0.00011542577452508812
    coeff6=-2.150293392583585e-7
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[8],memslots[2],memslots[3])
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-2.0909266321867902e-8
    coeff2=-2.7287925158899476e-7
    coeff3=-8.505783380079528e-8
    coeff4=-2.4699349369528174e-7
    coeff5=7.547096805172991e-6
    coeff6=-0.00011038627353830762
    coeff7=0.9999999715919681
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-0.00012498199105449398
    coeff2=0.9989757606293848
    coeff3=0.04754175140634288
    coeff4=0.002515599835204672
    coeff5=8.584302495182615e-5
    coeff6=5.619813578008144e-7
    coeff7=7.644984585045822e-9
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots[9],memslots[2],memslots[3])
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 3
    # Computing Ba8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.40308943704974315
    coeff2=-1.1533961194636121e-8
    coeff3=-5.03389309270026e-9
    coeff4=-1.635242146052482e-8
    coeff5=1.345288563528515e-7
    coeff6=8.355043856393377e-7
    coeff7=-0.00012492568018016112
    coeff8=0.9999999932294968
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=-1.0645153851732277e-9
    coeff2=0.9987013218201104
    coeff3=0.052949897775355415
    coeff4=0.003442617642989778
    coeff5=0.00011084212403663171
    coeff6=-3.8270257669926835e-7
    coeff7=-1.4340108650239077e-8
    coeff8=9.818180030024404e-11
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B8 with operation: mult
    mul!(memslots[10],memslots[2],memslots[3])
    # Deallocating Ba8 in slot 2
    # Deallocating Bb8 in slot 3
    # Computing T2k10 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=1.00000000175563
    coeff2=1.0000017846980391
    coeff3=0.49769531969924596
    coeff4=0.16600417761584907
    coeff5=0.03838601189631294
    coeff6=0.006185920919246356
    coeff7=0.0010600598008071952
    coeff8=8.103038123336939e-5
    coeff9=1.7925485378650857e-6
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9] .+ coeff9.*memslots[10]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    # Deallocating B7 in slot 9
    # Deallocating B8 in slot 10
    return memslots[1] # Returning T2k10
end

