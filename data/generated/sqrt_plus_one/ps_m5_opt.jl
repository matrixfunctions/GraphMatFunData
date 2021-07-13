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

@inline function ps_m5_opt(A)
    return ps_m5_opt!(copy(A))
end

@inline function ps_m5_opt!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=8
    memslots=Vector{Matrix{T}}(undef,max_memslots)
    n=size(A,1)
    for j=1:max_memslots
        memslots[j]=Matrix{T}(undef,n,n)
    end
    # The first slots are precomputed nodes [:A]
    memslots[1]=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 T2k8
    # Computing Bb2 = x*I+x*A
    coeff1=-0.15614734226906415
    coeff2=1.0138950449120168
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=-0.1561459961705319
    coeff2=1.0138962847063822
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=-1.1763384175098646
    coeff2=0.14819920360945626
    coeff3=0.7004279880959638
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=-0.15896658336808092
    coeff2=1.03474409787952
    coeff3=-0.6899181523788153
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots[5],memslots[2],memslots[3])
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.32868543911603565
    coeff2=0.04948779256327628
    coeff3=1.1645163087525163
    coeff4=0.33379574828605685
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=0.2715485387807492
    coeff2=0.4515173495475029
    coeff3=-0.45819098373171513
    coeff4=0.5179573824690393
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[2],memslots[3])
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.1784350608747074
    coeff2=-0.9109875948763699
    coeff3=0.20063142586799526
    coeff4=-1.6463200998393868
    coeff5=0.8066435555717062
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.9977699469177798
    coeff2=0.7154077419527756
    coeff3=-1.0589444472962166
    coeff4=0.32574547503223705
    coeff5=0.04488529740138668
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[2],memslots[3])
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.7570005286093737
    coeff2=0.007598464345057781
    coeff3=-0.8669783274817463
    coeff4=-0.14407160251709974
    coeff5=-0.6983909370796325
    coeff6=1.3011465317844473
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.6207474625776077
    coeff2=-0.8881294075403886
    coeff3=0.24425285144993336
    coeff4=-0.9312258295819543
    coeff5=-1.0814610150998485
    coeff6=1.0407649022964955
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[8],memslots[2],memslots[3])
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing T2k8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.2989395453295784
    coeff2=2.080346787066744
    coeff3=-1.5251262809722197
    coeff4=0.004244973897574398
    coeff5=-0.0326540011916379
    coeff6=0.8226732267898158
    coeff7=0.7171476120487754
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    return memslots[1] # Returning T2k8
end

