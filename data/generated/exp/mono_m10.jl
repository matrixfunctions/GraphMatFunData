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

@inline function mono_m10(A)
    return mono_m10!(copy(A))
end

@inline function mono_m10!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=11
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
    # Computation order: B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 T2k13
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing B3 with operation: mult
    mul!(memslots[3],memslots[2],memslots[1])
    # Computing B4 with operation: mult
    mul!(memslots[4],memslots[3],memslots[1])
    # Computing B5 with operation: mult
    mul!(memslots[5],memslots[4],memslots[1])
    # Computing B6 with operation: mult
    mul!(memslots[6],memslots[5],memslots[1])
    # Computing B7 with operation: mult
    mul!(memslots[7],memslots[6],memslots[1])
    # Computing B8 with operation: mult
    mul!(memslots[8],memslots[7],memslots[1])
    # Computing B9 with operation: mult
    mul!(memslots[9],memslots[8],memslots[1])
    # Computing B10 with operation: mult
    mul!(memslots[10],memslots[9],memslots[1])
    # Computing B11 with operation: mult
    mul!(memslots[11],memslots[10],memslots[1])
    # Computing T2k13 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8+x*B9+x*B10+x*B11
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=0.16666666666666666
    coeff5=0.041666666666666664
    coeff6=0.008333333333333333
    coeff7=0.001388888888888889
    coeff8=0.0001984126984126984
    coeff9=2.48015873015873e-5
    coeff10=2.7557319223985893e-6
    coeff11=2.755731922398589e-7
    coeff12=2.505210838544172e-8
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4] .+ coeff6.*memslots[5] .+ coeff7.*memslots[6] .+ coeff8.*memslots[7] .+ coeff9.*memslots[8] .+ coeff10.*memslots[9] .+ coeff11.*memslots[10] .+ coeff12.*memslots[11]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 4
    # Deallocating B5 in slot 5
    # Deallocating B6 in slot 6
    # Deallocating B7 in slot 7
    # Deallocating B8 in slot 8
    # Deallocating B9 in slot 9
    # Deallocating B10 in slot 10
    # Deallocating B11 in slot 11
    return memslots[1] # Returning T2k13
end

