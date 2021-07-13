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


# Compute X <- a X + b Y.
function matfun_axpby!(X,a,b,Y)
    m,n=size(X)
    if ~(a isa ValueOne)
        rmul!(X,a)
    end
    @inbounds for i=1:m
        @inbounds for j=1:n
            if (b isa ValueOne)
                X[i,j]+=Y[i,j]
            else
                X[i,j]+=b*Y[i,j]
            end
        end
    end
end

@inline function exp_sastre_m6(A)
    return exp_sastre_m6!(copy(A))
end

@inline function exp_sastre_m6!(A)
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
    # Computation order: B2 Bb3 B3 Ba4 Bb4 B4 Ba7_5 Bb5_3 B5 Ba6 Bb6 B6 Bb7 B7 T2k9
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing Bb3 = x*A+x*B2
    coeff1=0.001052151783051235
    coeff2=0.0004675683454147702
    memslots[3] .= coeff1.*memslots[1] .+ coeff2.*memslots[2]
    # Computing B3 with operation: mult
    mul!(memslots[4],memslots[2],memslots[3])
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*A+x*B2+x*B3
    coeff1=0.2868706220817633
    coeff2=-0.03289442879547955
    coeff3=1.0
    memslots[3] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[4]
    # Computing Bb4 = x*B2+x*B3
    coeff1=0.05317514832355802
    coeff2=1.0
    memslots[5] .= coeff1.*memslots[2] .+ coeff2.*memslots[4]
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[3],memslots[5])
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 5
    # Computing Ba7_5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=3.09646797193604
    coeff2=0.772360321294401
    coeff3=0.1673139636901279
    coeff4=7.922322450524197
    coeff5=1.0
    # Smart lincomb recycle B3
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[4] .+ coeff5.*memslots[6]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Deallocating B4 in slot 6
    # Computing Bb5_3 = x*A+x*B2
    coeff1=0.002688394980266927
    coeff2=0.0004675683454147702
    memslots[3] .= coeff1.*memslots[1] .+ coeff2.*memslots[2]
    # Computing B5 with operation: mult
    mul!(memslots[5],memslots[2],memslots[3])
    # Deallocating Bb5_3 in slot 3
    # Computing Ba6 = x*A+x*B2+x*B5
    coeff1=0.39689859154115
    coeff2=0.02219811707032801
    coeff3=1.0
    memslots[3] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[5]
    # Computing Bb6 = x*B2+x*B5
    coeff1=0.0277140002806296
    coeff2=1.0
    memslots[6] .= coeff1.*memslots[2] .+ coeff2.*memslots[5]
    # Computing B6 with operation: mult
    mul!(memslots[7],memslots[3],memslots[6])
    # Deallocating Ba6 in slot 3
    # Deallocating Bb6 in slot 6
    # Computing Bb7 = x*A+x*B2+x*B5+x*B6
    coeff1=0.3229486011362678
    coeff2=0.08092036376147299
    coeff3=1.930814505527068
    coeff4=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[5] .+ coeff4.*memslots[7]
    # Deallocating B2 in slot 2
    # Deallocating B5 in slot 5
    # Deallocating B6 in slot 7
    # Computing B7 with operation: mult
    mul!(memslots[2],memslots[4],memslots[1])
    # Deallocating Ba7_5 in slot 4
    # Deallocating Bb7 in slot 1
    # Computing T2k9 = x*I+x*B7
    coeff1=1.0
    coeff2=1.0
    # Smart lincomb recycle B7
    memslots[2] .= coeff2.*memslots[2]
    mul!(memslots[2],true,I*coeff1,true,true)
    return memslots[2] # Returning T2k9
end

