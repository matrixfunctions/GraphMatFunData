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

@inline function exp_ps_m7_opt_rho6_4(A)
    return exp_ps_m7_opt_rho6_4!(copy(A))
end

@inline function exp_ps_m7_opt_rho6_4!(A)
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
    coeff1=0.42031771622120523
    coeff2=0.309252105381124
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=0.42031771622120523
    coeff2=0.309252105381124
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=0.09654809779701808
    coeff2=0.9314722195252554
    coeff3=0.19731022777550086
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=0.9062397883450928
    coeff2=0.11339252602614255
    coeff3=0.12302744324936216
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots[5],memslots[2],memslots[3])
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.06178357797354583
    coeff2=0.13477081084608053
    coeff3=0.7329752440067754
    coeff4=0.5453304515439815
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=0.7648819021317048
    coeff2=0.3110254108862011
    coeff3=0.13872927228389606
    coeff4=0.031072070853003377
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[2],memslots[3])
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.0037226968956481066
    coeff2=0.3027942237668864
    coeff3=0.16951280499357504
    coeff4=0.7173684223758208
    coeff5=0.6151562256858113
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.4984122280039027
    coeff2=0.6168934079023138
    coeff3=0.32464232547450916
    coeff4=0.11022586483948126
    coeff5=0.00562561032811143
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[2],memslots[3])
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-2.47688825760804e-9
    coeff2=0.0004847383569867066
    coeff3=-0.008759615364757254
    coeff4=-0.024971360914667298
    coeff5=0.20141416233600964
    coeff6=0.9838034227912541
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-7.332920860083383e-8
    coeff2=7.255106883862399e-10
    coeff3=1.599720660809858e-11
    coeff4=5.386947932212966e-11
    coeff5=3.696432230206371e-12
    coeff6=1.7604767541951467e-14
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[8],memslots[2],memslots[3])
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=5.241383560289876e-5
    coeff2=0.014487073789250026
    coeff3=0.01054251273835961
    coeff4=-0.12590635459023522
    coeff5=-0.35076815913299725
    coeff6=0.9784008568588043
    coeff7=4.756683002701541e-7
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-0.004732216406789163
    coeff2=2.4676441956771072e-5
    coeff3=1.805471807800756e-5
    coeff4=6.965400723989983e-6
    coeff5=1.4212591399083807e-6
    coeff6=1.2164967481538e-7
    coeff7=1.0000000002008125
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots[9],memslots[2],memslots[3])
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 3
    # Computing Ba8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.9900217594371604
    coeff2=-0.0426544708742059
    coeff3=0.0895038182070465
    coeff4=0.4332397850086141
    coeff5=0.6230790774477177
    coeff6=0.836798518267107
    coeff7=6.715240076479665e-7
    coeff8=0.03410244422083362
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.09171588796170577
    coeff2=0.01956469387438758
    coeff3=0.01590736865892108
    coeff4=0.010821193755227018
    coeff5=0.002975380215829149
    coeff6=0.005577510585906947
    coeff7=-7.0435154500906185e-6
    coeff8=1.0199299719677057
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B8 with operation: mult
    mul!(memslots[10],memslots[2],memslots[3])
    # Deallocating Ba8 in slot 2
    # Deallocating Bb8 in slot 3
    # Computing T2k10 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=0.7602626115113443
    coeff2=0.2868135868803526
    coeff3=0.36416792067469345
    coeff4=0.1470803112875071
    coeff5=0.301249740245868
    coeff6=0.05977870243645669
    coeff7=3.4088571170380866e-7
    coeff8=0.00148688702927498
    coeff9=1.021358342405418
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

