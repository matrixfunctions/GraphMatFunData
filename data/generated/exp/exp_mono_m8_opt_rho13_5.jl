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

@inline function exp_mono_m8_opt_rho13_5(A)
    return exp_mono_m8_opt_rho13_5!(copy(A))
end

@inline function exp_mono_m8_opt_rho13_5!(A)
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
    # Computation order: B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 Ba7 Bb7 B7 Ba8 Bb8 B8 Ba9 Bb9 B9 T2k11
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=-0.003584018414055765
    coeff2=0.16376260367565243
    coeff3=0.04917059370539749
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[2]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=0.30528232877254763
    coeff2=0.23746768414771152
    coeff3=0.0023273769915510135
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots[5],memslots[3],memslots[4])
    # Deallocating Ba3 in slot 3
    # Deallocating Bb3 in slot 4
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=0.0005515398938947843
    coeff2=0.02533967378579724
    coeff3=0.017067262331527198
    coeff4=0.9638185168573307
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.008872840493334803
    coeff2=0.25571769129339744
    coeff3=0.0037256194747560907
    coeff4=0.002052240135641998
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[3],memslots[4])
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 4
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.0012630062267903097
    coeff2=-0.0010962696861478022
    coeff3=0.00022839563887310042
    coeff4=0.013673005502563394
    coeff5=0.8216219421030205
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.8293727566445452
    coeff2=0.1675033881979259
    coeff3=0.0030839993650848517
    coeff4=0.0031897875289976277
    coeff5=3.234759791150232e-5
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[3],memslots[4])
    # Deallocating Ba5 in slot 3
    # Deallocating Bb5 in slot 4
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=4.437908455542293e-5
    coeff2=0.0003258145596094315
    coeff3=0.00022081927999800417
    coeff4=-0.13039438638797307
    coeff5=0.8475672106088131
    coeff6=0.8145974867536132
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-0.10275230374603168
    coeff2=0.25648676493721656
    coeff3=0.005628536962635935
    coeff4=0.008329230127484285
    coeff5=0.00021583757875408033
    coeff6=-1.0968604718569746e-7
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[8],memslots[3],memslots[4])
    # Deallocating Ba6 in slot 3
    # Deallocating Bb6 in slot 4
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-2.8201037195155804e-6
    coeff2=4.8728643050646066e-6
    coeff3=2.2723646529748017e-6
    coeff4=-0.0061294802470251075
    coeff5=0.0688513777099278
    coeff6=-0.12038470670165373
    coeff7=0.9948034699390058
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.004650329563627054
    coeff2=0.2598813154789629
    coeff3=0.006290014251994485
    coeff4=0.010576806709956765
    coeff5=0.0003949349245182783
    coeff6=8.856309878298438e-6
    coeff7=7.553245375703866e-9
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots[9],memslots[3],memslots[4])
    # Deallocating Ba7 in slot 3
    # Deallocating Bb7 in slot 4
    # Computing Ba8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=1.0002604649982476
    coeff2=0.2422828857323093
    coeff3=0.022001533431922154
    coeff4=0.15933029209747282
    coeff5=0.03184541272812379
    coeff6=0.007635550247315368
    coeff7=0.0008615033049494069
    coeff8=3.134369854398072e-5
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing Bb8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=1.0002604649982476
    coeff2=0.2422828857323093
    coeff3=0.022001533431922154
    coeff4=0.15933029209747282
    coeff5=0.03184541272812379
    coeff6=0.007635550247315368
    coeff7=0.0008615033049494069
    coeff8=3.134369854398072e-5
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing B8 with operation: mult
    mul!(memslots[10],memslots[3],memslots[4])
    # Deallocating Ba8 in slot 3
    # Deallocating Bb8 in slot 4
    # Computing Ba9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=0.1346603332119568
    coeff2=-5.9198518169123926e-5
    coeff3=-8.660244928553117e-6
    coeff4=-8.731780915688725e-5
    coeff5=-2.1711606510903095e-5
    coeff6=-5.9723490061415135e-6
    coeff7=-7.70252692547925e-7
    coeff8=-3.1855167508881976e-8
    coeff9=1.0081137407871135
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9] .+ coeff9.*memslots[10]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing Bb9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=-0.06540541282818779
    coeff2=6.846710342127941e-5
    coeff3=2.738382227960167e-6
    coeff4=-6.3730674335546956e-6
    coeff5=-5.739005263969626e-6
    coeff6=-2.157568588204766e-6
    coeff7=-3.3636703984819904e-7
    coeff8=-1.5592924041261544e-8
    coeff9=0.9946441317159572
    # Smart lincomb recycle B7
    memslots[9] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9] .+ coeff9.*memslots[10]
    mul!(memslots[9],true,I*coeff1,true,true)
    # Computing B9 with operation: mult
    mul!(memslots[4],memslots[3],memslots[9])
    # Deallocating Ba9 in slot 3
    # Deallocating Bb9 in slot 9
    # Computing T2k11 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B8+x*B9
    coeff1=0.008776543046053789
    coeff2=-1.3078408314728238e-5
    coeff3=-9.364829337290504e-7
    coeff4=-4.881187067247219e-6
    coeff5=-6.561388347938557e-7
    coeff6=-1.0279110807751132e-7
    coeff7=-5.458270135394687e-9
    coeff8=-0.06726697857924235
    coeff9=0.9964850858948753
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[10] .+ coeff9.*memslots[4]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    # Deallocating B8 in slot 10
    # Deallocating B9 in slot 4
    return memslots[1] # Returning T2k11
end

