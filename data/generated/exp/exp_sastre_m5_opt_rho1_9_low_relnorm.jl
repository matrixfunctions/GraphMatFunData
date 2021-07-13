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

@inline function exp_sastre_m5_opt_rho1_9_low_relnorm(A)
    return exp_sastre_m5_opt_rho1_9_low_relnorm!(copy(A))
end

@inline function exp_sastre_m5_opt_rho1_9_low_relnorm!(A)
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
    coeff1=-0.31551655352795516
    coeff2=0.4801289831900131
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=-0.31551655352795516
    coeff2=0.4801289831900131
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=-0.357355571021064
    coeff2=0.4105933675776958
    coeff3=0.09852026651116912
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=0.053558315750694814
    coeff2=-0.22014016998509814
    coeff3=0.922091412659715
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots[5],memslots[2],memslots[3])
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.004620908141128493
    coeff2=0.055057556018866645
    coeff3=-0.351290300034291
    coeff4=0.9607057168526784
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=0.011326704451019963
    coeff2=0.002770498744594938
    coeff3=0.00045550004840354903
    coeff4=0.00012416495659100344
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[2],memslots[3])
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.11176905543585865
    coeff2=0.5925762886270637
    coeff3=0.4621864853147797
    coeff4=0.038392790243209074
    coeff5=0.24031719666304735
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.1215251552594178
    coeff2=0.2145595643740629
    coeff3=0.04592597180399795
    coeff4=0.021375614386974554
    coeff5=0.512024337751819
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[2],memslots[3])
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=1.11114318963272
    coeff2=0.4876551310270585
    coeff3=0.3710914185617396
    coeff4=0.08701143872628835
    coeff5=5.180154142825969
    coeff6=0.6495647957423022
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.9070495211990466
    coeff2=0.5443460967617421
    coeff3=0.5323384854781406
    coeff4=0.12934538808055288
    coeff5=5.16322755510199
    coeff6=0.6107117158792567
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[8],memslots[2],memslots[3])
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing T2k8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-0.0818796601588853
    coeff2=0.003683007236493971
    coeff3=0.008733812199626457
    coeff4=0.005536448121469381
    coeff5=0.003015592881979151
    coeff6=0.01836716605825683
    coeff7=0.9780133820655359
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

