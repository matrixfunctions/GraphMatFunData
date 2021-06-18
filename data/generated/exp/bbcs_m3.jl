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

@inline function bbcs_m3(A)
    return bbcs_m3!(copy(A))
end

@inline function bbcs_m3!(A)
    T=promote_type(eltype(A),ComplexF64) # Make it work for many 'bigger' types (matrices and scalars)
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
    # Computation order: B2 Bb3 B3 Ba4 Bb4 B4 T2k6
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing Bb3 = x*A+x*B2
    coeff1=0.10775 + 0.0im
    coeff2=-0.0 - 0.026939068735988708im
    memslots[3] .= coeff1.*memslots[1] .+ coeff2.*memslots[2]
    # Computing B3 with operation: mult
    mul!(memslots[4],memslots[2],memslots[3])
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*B2+x*B3
    coeff1=0.0 + 0.6632100444166243im
    coeff2=1.0 + 0.0im
    memslots[3] .= coeff1.*memslots[2] .+ coeff2.*memslots[4]
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=0.0 + 0.5496085391143601im
    coeff2=0.1620095284677366 + 0.0im
    coeff3=0.0 - 0.014179818052118045im
    coeff4=-0.034159539168921116 + 0.0im
    # Smart lincomb recycle B3
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[4]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[5],memslots[3],memslots[4])
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 4
    # Computing T2k6 = x*I+x*A+x*B2+x*B4
    coeff1=1.0 + 0.0im
    coeff2=0.0 - 0.9999999999999923im
    coeff3=-0.13549409636220702 + 0.0im
    coeff4=1.0 + 0.0im
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[5]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B4 in slot 5
    return memslots[1] # Returning T2k6
end
