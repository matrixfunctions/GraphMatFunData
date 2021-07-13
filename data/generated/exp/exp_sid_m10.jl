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

@inline function exp_sid_m10(A)
    return exp_sid_m10!(copy(A))
end

@inline function exp_sid_m10!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=9
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
    # Computation order: Bb2 Ba2 B2 Ba3 B3 Ba4 B4 Ba5 B5 Bb6 B6 Ba7 Bb7 B7 Ba8 Bb8_7 B8 Ba9 Bb9 B9 B10 B11 T2k13
    # Computing Bb2 = x*I+x*A
    coeff1=0.0
    coeff2=0.125
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=0.0
    coeff2=0.125
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*A+x*B2
    coeff1=0.125
    coeff2=0.0
    memslots[2] .= coeff1.*memslots[1] .+ coeff2.*memslots[4]
    # Computing B3 with operation: mult
    mul!(memslots[3],memslots[2],memslots[4])
    # Deallocating Ba3 in slot 2
    # Computing Ba4 = x*A+x*B3
    coeff1=0.125
    coeff2=0.0
    memslots[2] .= coeff1.*memslots[1] .+ coeff2.*memslots[3]
    # Computing B4 with operation: mult
    mul!(memslots[5],memslots[2],memslots[3])
    # Deallocating Ba4 in slot 2
    # Computing Ba5 = x*A+x*B4
    coeff1=0.125
    coeff2=0.0
    memslots[2] .= coeff1.*memslots[1] .+ coeff2.*memslots[5]
    # Computing B5 with operation: mult
    mul!(memslots[6],memslots[2],memslots[5])
    # Deallocating Ba5 in slot 2
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=1.6421017771235788e-7
    coeff2=6.204734935438909e-8
    coeff3=2.957106114715868e-9
    coeff4=1.556371639324141e-10
    coeff5=1.556371639324141e-11
    memslots[2] .= coeff1.*memslots[1] .+ coeff2.*memslots[4] .+ coeff3.*memslots[3] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    # Computing B6 with operation: mult
    mul!(memslots[7],memslots[6],memslots[2])
    # Deallocating Bb6 in slot 2
    # Computing Ba7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.7439481579382581
    coeff2=0.4155284057336423
    coeff3=0.02479095151834799
    coeff4=0.001283057135586989
    coeff5=3.501669195497238e-5
    coeff6=1.0
    memslots[2] .= coeff1.*memslots[1] .+ coeff2.*memslots[4] .+ coeff3.*memslots[3] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    # Computing Bb7 = x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.03306559506631931
    coeff2=0.002630043177655382
    coeff3=0.0002100333647757715
    coeff4=3.7537107416419e-5
    coeff5=1.0
    memslots[8] .= coeff1.*memslots[4] .+ coeff2.*memslots[3] .+ coeff3.*memslots[5] .+ coeff4.*memslots[6] .+ coeff5.*memslots[7]
    # Computing B7 with operation: mult
    mul!(memslots[9],memslots[2],memslots[8])
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 8
    # Computing Ba8 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=1.3883617476073524
    coeff2=2.991654767354374
    coeff3=0.2857950268422422
    coeff4=0.03005135891320298
    coeff5=0.002742336655922557
    coeff6=61.75954247606858
    coeff7=1.0
    memslots[2] .= coeff1.*memslots[1] .+ coeff2.*memslots[4] .+ coeff3.*memslots[3] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[9]
    # Computing Bb8_7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-3.2977952779222e-5
    coeff2=0.008139086096860678
    coeff3=0.001121744731945438
    coeff4=9.027588625491207e-5
    coeff5=8.572383602707347e-6
    coeff6=1.0
    # Smart lincomb recycle B6
    memslots[7] .= coeff1.*memslots[1] .+ coeff2.*memslots[4] .+ coeff3.*memslots[3] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    # Computing B8 with operation: mult
    mul!(memslots[8],memslots[2],memslots[7])
    # Deallocating Ba8 in slot 2
    # Deallocating Bb8_7 in slot 7
    # Computing Ba9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B8
    coeff1=1.0
    coeff2=0.125
    coeff3=0.5029302610017967
    coeff4=0.07705596948494946
    coeff5=0.004985549176118462
    coeff6=6.263526066651383e-5
    coeff7=1.0
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[3] .+ coeff5.*memslots[5] .+ coeff6.*memslots[6] .+ coeff7.*memslots[8]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B8
    coeff1=1.0
    coeff2=0.125
    coeff3=0.5029302610017967
    coeff4=0.07705596948494946
    coeff5=0.004985549176118462
    coeff6=6.263526066651383e-5
    coeff7=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[3] .+ coeff5.*memslots[5] .+ coeff6.*memslots[6] .+ coeff7.*memslots[8]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 5
    # Deallocating B5 in slot 6
    # Deallocating B8 in slot 8
    # Computing B9 with operation: mult
    mul!(memslots[3],memslots[2],memslots[1])
    # Deallocating Ba9 in slot 2
    # Deallocating Bb9 in slot 1
    # Computing B10 with operation: mult
    mul!(memslots[1],memslots[3],memslots[3])
    # Deallocating B9 in slot 3
    # Computing B11 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Deallocating B10 in slot 1
    # Computing T2k13 = x*B7+x*B11
    coeff1=0.0
    coeff2=1.0
    # Smart lincomb recycle B7
    memslots[9] .= coeff1.*memslots[9] .+ coeff2.*memslots[2]
    # Deallocating B11 in slot 2
    return memslots[9] # Returning T2k13
end

