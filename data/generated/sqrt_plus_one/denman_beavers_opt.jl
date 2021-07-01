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

@inline function denman_beavers_opt(A)
    return denman_beavers_opt!(copy(A))
end

@inline function denman_beavers_opt!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=7
    memslots=Vector{Matrix{T}}(undef,max_memslots)
    n=size(A,1)
    for j=1:max_memslots
        memslots[j]=Matrix{T}(undef,n,n)
    end
    # The first slots are precomputed nodes [:A]
    memslots[1]=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 T2k7
    # Computing Bb2 = x*I+x*A
    coeff1=1.0
    coeff2=-0.011141376535527443
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=1.0
    coeff2=0.7610413081657074
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: ldiv
    memslots[4]=memslots[3]\memslots[2]
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=1.0
    coeff2=-0.07127247848266706
    coeff3=0.9541397430184381
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=1.0
    coeff2=-0.02575361014723645
    coeff3=-0.004282806676350282
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B3 with operation: ldiv
    memslots[5]=memslots[2]\memslots[3]
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=1.0
    coeff2=0.5152503668289159
    coeff3=-0.10150239642416546
    coeff4=0.016950835467818528
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=1.0
    coeff2=0.03845623901890617
    coeff3=0.026744368550352622
    coeff4=-0.029043605652931636
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B4 with operation: ldiv
    memslots[6]=memslots[2]\memslots[3]
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.49421282325538585
    coeff2=0.05943417585870563
    coeff3=0.5017024742690791
    coeff4=0.5276435956563207
    coeff5=1.0
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=1.0
    coeff2=0.12762852235945304
    coeff3=-0.0015101375211018255
    coeff4=-0.14486471812413823
    coeff5=0.01583234295751802
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B5 with operation: ldiv
    memslots[7]=memslots[2]\memslots[3]
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing T2k7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.28457285753903816
    coeff2=0.11495283295141033
    coeff3=-0.0344930890633602
    coeff4=0.34606871085737384
    coeff5=-0.006442438999085179
    coeff6=1.4670072325444539
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    return memslots[1] # Returning T2k7
end

