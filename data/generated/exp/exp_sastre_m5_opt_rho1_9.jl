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

@inline function exp_sastre_m5_opt_rho1_9(A)
    return exp_sastre_m5_opt_rho1_9!(copy(A))
end

@inline function exp_sastre_m5_opt_rho1_9!(A)
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
    coeff1=-0.29597130380221776
    coeff2=0.4830650480818557
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=-0.29597130380221776
    coeff2=0.4830650480818557
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=-0.35706341164511884
    coeff2=0.41479581112660663
    coeff3=0.08232725511150228
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=0.052794118729278604
    coeff2=-0.22109848126655623
    coeff3=0.9222178419369582
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots[5],memslots[2],memslots[3])
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.005411040008503428
    coeff2=0.0479465268958071
    coeff3=-0.3539095385982056
    coeff4=0.9600894604220988
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=0.014197999026813905
    coeff2=0.0023989722622025037
    coeff3=0.0004504909091063826
    coeff4=0.00011985947592822976
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[2],memslots[3])
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.10677883331398988
    coeff2=0.6584173995757578
    coeff3=0.4164266931927851
    coeff4=0.03716863529834358
    coeff5=0.2855172295604338
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.12273251114319797
    coeff2=0.1799587581454721
    coeff3=0.043513380168846295
    coeff4=0.02122483727644592
    coeff5=0.5818565692388099
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[2],memslots[3])
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=1.0985965640823765
    coeff2=0.45359605154012905
    coeff3=0.3807951389639113
    coeff4=0.09228091229564235
    coeff5=5.179587866680269
    coeff6=0.6801751580137498
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.9068877097254958
    coeff2=0.5793395460985866
    coeff3=0.5479499415010913
    coeff4=0.13249960107719003
    coeff5=5.153379110061045
    coeff6=0.6370429665254512
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[8],memslots[2],memslots[3])
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing T2k8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-0.050949343525589726
    coeff2=-0.019581395624785863
    coeff3=-0.016433672152666983
    coeff4=-0.0018851674132650941
    coeff5=0.0009698855603010127
    coeff6=0.014316598391468338
    coeff7=0.9728079294146128
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

