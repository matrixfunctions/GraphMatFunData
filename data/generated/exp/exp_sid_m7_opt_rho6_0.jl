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

@inline function exp_sid_m7_opt_rho6_0(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_sid_m7_opt_rho6_0!(A_copy)
end

@inline function exp_sid_m7_opt_rho6_0!(A)
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
    coeff1=-0.0674550728089079
    coeff2=1.1639373098746155
    memslots2 .= coeff2.*memslots1
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=-0.0674550728089079
    coeff2=1.1639373098746155
    memslots3 .= coeff2.*memslots1
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=-0.05790075880684695
    coeff2=1.050095804216935
    coeff3=0.02211891467634549
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=-0.000951838725229865
    coeff2=-0.035990141876148504
    coeff3=1.0505773512582828
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots5,memslots2,memslots3)
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.022273287604033215
    coeff2=0.999419730173801
    coeff3=0.04428302676636615
    coeff4=0.0013986894432222274
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.0003564977821414783
    coeff2=0.006790605847297012
    coeff3=-0.0863934611339156
    coeff4=0.9986696579293114
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots6,memslots2,memslots3)
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.00992506274997938
    coeff2=1.000247094423774
    coeff3=0.010965295621586402
    coeff4=-0.0007441489467712597
    coeff5=3.3836127966022704e-6
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.00010751494099623147
    coeff2=0.00045481726705312207
    coeff3=0.0010096879534240845
    coeff4=-0.027269903163941744
    coeff5=0.999827547479686
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots7,memslots2,memslots3)
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=2.8512702775691487e-6
    coeff2=8.323858094672198e-6
    coeff3=1.3339497680252084e-5
    coeff4=-0.0005347112859929923
    coeff5=-0.009928519927039425
    coeff6=1.0000008787568422
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-3.9152959471223135e-5
    coeff2=8.334641611255522e-8
    coeff3=7.3222974116279275e-9
    coeff4=1.8759928836648602e-10
    coeff5=1.7224589038431254e-12
    coeff6=1.6134727895776298e-13
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots8,memslots2,memslots3)
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-0.0006322010382109642
    coeff2=5.966629547283877
    coeff3=0.3119218262051087
    coeff4=0.01788182758740769
    coeff5=0.0005913130477887621
    coeff6=4.309674736291994e-5
    coeff7=1.1067710016183396
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-0.034743137725474355
    coeff2=0.12434316511921449
    coeff3=0.04622622385713797
    coeff4=0.006600522551149245
    coeff5=0.000256268554613047
    coeff6=4.892795851939188e-5
    coeff7=1.1438201195458555
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots9,memslots2,memslots3)
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 3
    # Computing Ba8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.20063211995463168
    coeff2=11.100148209740121
    coeff3=3.0262688473369668
    coeff4=0.39898152652406355
    coeff5=0.01387684466194758
    coeff6=0.002742528074882574
    coeff7=61.75609557189336
    coeff8=1.1419333121703887
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.013739397806079797
    coeff2=-0.034167351013704314
    coeff3=-0.00017439769298246415
    coeff4=0.0006895671006188639
    coeff5=3.0262187841524036e-5
    coeff6=3.93368482319488e-5
    coeff7=0.9924325575605532
    coeff8=-7.082137587033525e-5
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B8 with operation: mult
    mul!(memslots10,memslots2,memslots3)
    # Deallocating Ba8 in slot 2
    # Deallocating Bb8 in slot 3
    # Computing T2k10 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=0.9943018399475028
    coeff2=0.9623400714214185
    coeff3=0.621285209008396
    coeff4=0.22263803654723124
    coeff5=0.04001562250761491
    coeff6=0.0017973977825017301
    coeff7=0.0002723523751212443
    coeff8=0.012856181363394585
    coeff9=0.9918580039315373
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

