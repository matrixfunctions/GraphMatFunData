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

@inline function ps(A)
    return ps!(copy(A))
end

@inline function ps!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=5
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
    # Computation order: B2 B3 Bb4 B4 Bb5 B5 T2k7
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing B3 with operation: mult
    mul!(memslots[3],memslots[2],memslots[1])
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.0205078125
    coeff2=0.01611328125
    coeff3=-0.013092041015625
    coeff4=0.0109100341796875
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[5],memslots[3],memslots[4])
    # Deallocating Bb4 in slot 4
    # Computing Bb5 = x*I+x*A+x*B2+x*B4
    coeff1=0.0625
    coeff2=-0.0390625
    coeff3=0.02734375
    coeff4=1.0
    # Smart lincomb recycle B4
    memslots[5] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5]
    mul!(memslots[5],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[4],memslots[3],memslots[5])
    # Deallocating B3 in slot 3
    # Deallocating Bb5 in slot 5
    # Computing T2k7 = x*I+x*A+x*B2+x*B5
    coeff1=1.0
    coeff2=0.5
    coeff3=-0.125
    coeff4=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[4]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B5 in slot 4
    return memslots[1] # Returning T2k7
end

