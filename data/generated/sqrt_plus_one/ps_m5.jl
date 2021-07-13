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

@inline function ps_m5(A)
    return ps_m5!(copy(A))
end

@inline function ps_m5!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=6
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
    # Computation order: B2 B3 B4 Bb5 B5 Bb6 B6 T2k8
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing B3 with operation: mult
    mul!(memslots[3],memslots[2],memslots[1])
    # Computing B4 with operation: mult
    mul!(memslots[4],memslots[3],memslots[1])
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.013092041015625
    coeff2=0.0109100341796875
    coeff3=-0.009273529052734375
    coeff4=0.008008956909179688
    coeff5=-0.0070078372955322266
    memslots[5] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4]
    mul!(memslots[5],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[6],memslots[4],memslots[5])
    # Deallocating Bb5 in slot 5
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B5
    coeff1=-0.0390625
    coeff2=0.02734375
    coeff3=-0.0205078125
    coeff4=0.01611328125
    coeff5=1.0
    # Smart lincomb recycle B5
    memslots[6] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[6]
    mul!(memslots[6],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[5],memslots[4],memslots[6])
    # Deallocating B4 in slot 4
    # Deallocating Bb6 in slot 6
    # Computing T2k8 = x*I+x*A+x*B2+x*B3+x*B6
    coeff1=1.0
    coeff2=0.5
    coeff3=-0.125
    coeff4=0.0625
    coeff5=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[5]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B6 in slot 5
    return memslots[1] # Returning T2k8
end

