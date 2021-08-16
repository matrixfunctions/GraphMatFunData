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

@inline function exp_sid_m11(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_sid_m11!(A_copy)
end

@inline function exp_sid_m11!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=9
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
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: Bb2 Ba2 B2 Ba3 B3 Ba4 B4 Ba5 B5 Bb6 B6 Ba7 Bb7 B7 Ba8 Bb8_7 B8 Ba9 Bb9 B9 B10 B11 B12 T2k14
    # Computing Bb2 = x*I+x*A
    coeff1=0.0
    coeff2=0.0625
    memslots2 .= coeff2.*memslots1
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=0.0
    coeff2=0.0625
    memslots3 .= coeff2.*memslots1
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*A+x*B2
    coeff1=0.0625
    coeff2=0.0
    memslots2 .= coeff1.*memslots1 .+ coeff2.*memslots4
    # Computing B3 with operation: mult
    mul!(memslots3,memslots2,memslots4)
    # Deallocating Ba3 in slot 2
    # Computing Ba4 = x*A+x*B3
    coeff1=0.0625
    coeff2=0.0
    memslots2 .= coeff1.*memslots1 .+ coeff2.*memslots3
    # Computing B4 with operation: mult
    mul!(memslots5,memslots2,memslots3)
    # Deallocating Ba4 in slot 2
    # Computing Ba5 = x*A+x*B4
    coeff1=0.0625
    coeff2=0.0
    memslots2 .= coeff1.*memslots1 .+ coeff2.*memslots5
    # Computing B5 with operation: mult
    mul!(memslots6,memslots2,memslots5)
    # Deallocating Ba5 in slot 2
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=8.210508885617894e-8
    coeff2=6.204734935438909e-8
    coeff3=2.957106114715868e-9
    coeff4=1.556371639324141e-10
    coeff5=1.556371639324141e-11
    memslots2 .= coeff1.*memslots1 .+ coeff2.*memslots4 .+ coeff3.*memslots3 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    # Computing B6 with operation: mult
    mul!(memslots7,memslots6,memslots2)
    # Deallocating Bb6 in slot 2
    # Computing Ba7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.37197407896912904
    coeff2=0.4155284057336423
    coeff3=0.02479095151834799
    coeff4=0.001283057135586989
    coeff5=3.501669195497238e-5
    coeff6=1.0
    memslots2 .= coeff1.*memslots1 .+ coeff2.*memslots4 .+ coeff3.*memslots3 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    # Computing Bb7 = x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.03306559506631931
    coeff2=0.002630043177655382
    coeff3=0.0002100333647757715
    coeff4=3.7537107416419e-5
    coeff5=1.0
    memslots8 .= coeff1.*memslots4 .+ coeff2.*memslots3 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7
    # Computing B7 with operation: mult
    mul!(memslots9,memslots2,memslots8)
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 8
    # Computing Ba8 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.6941808738036762
    coeff2=2.991654767354374
    coeff3=0.2857950268422422
    coeff4=0.03005135891320298
    coeff5=0.002742336655922557
    coeff6=61.75954247606858
    coeff7=1.0
    memslots2 .= coeff1.*memslots1 .+ coeff2.*memslots4 .+ coeff3.*memslots3 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots9
    # Computing Bb8_7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-1.6488976389611e-5
    coeff2=0.008139086096860678
    coeff3=0.001121744731945438
    coeff4=9.027588625491207e-5
    coeff5=8.572383602707347e-6
    coeff6=1.0
    # Smart lincomb recycle B6
    memslots7 .= coeff1.*memslots1 .+ coeff2.*memslots4 .+ coeff3.*memslots3 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    # Computing B8 with operation: mult
    mul!(memslots8,memslots2,memslots7)
    # Deallocating Ba8 in slot 2
    # Deallocating Bb8_7 in slot 7
    # Computing Ba9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B8
    coeff1=1.0
    coeff2=0.0625
    coeff3=0.5029302610017967
    coeff4=0.07705596948494946
    coeff5=0.004985549176118462
    coeff6=6.263526066651383e-5
    coeff7=1.0
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots3 .+ coeff5.*memslots5 .+ coeff6.*memslots6 .+ coeff7.*memslots8
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B8
    coeff1=1.0
    coeff2=0.0625
    coeff3=0.5029302610017967
    coeff4=0.07705596948494946
    coeff5=0.004985549176118462
    coeff6=6.263526066651383e-5
    coeff7=1.0
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots3 .+ coeff5.*memslots5 .+ coeff6.*memslots6 .+ coeff7.*memslots8
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 5
    # Deallocating B5 in slot 6
    # Deallocating B8 in slot 8
    # Computing B9 with operation: mult
    mul!(memslots3,memslots2,memslots1)
    # Deallocating Ba9 in slot 2
    # Deallocating Bb9 in slot 1
    # Computing B10 with operation: mult
    mul!(memslots1,memslots3,memslots3)
    # Deallocating B9 in slot 3
    # Computing B11 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Deallocating B10 in slot 1
    # Computing B12 with operation: mult
    mul!(memslots1,memslots2,memslots2)
    # Deallocating B11 in slot 2
    # Computing T2k14 = x*B7+x*B12
    coeff1=0.0
    coeff2=1.0
    # Smart lincomb recycle B7
    memslots9 .= coeff1.*memslots9 .+ coeff2.*memslots1
    # Deallocating B12 in slot 1
    return memslots9 # Returning T2k14
end

