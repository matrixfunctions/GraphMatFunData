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

@inline function bbcs_m4(A)
    return bbcs_m4!(copy(A))
end

@inline function bbcs_m4!(A)
    T=promote_type(eltype(A),ComplexF64) # Make it work for many 'bigger' types (matrices and scalars)
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
    # Computation order: B2 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 T2k7
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing B3 with operation: mult
    mul!(memslots[3],memslots[2],memslots[1])
    # Computing Ba4 = x*A+x*B2+x*B3
    coeff1=0.0 + 0.13340427306445612im
    coeff2=0.020226020298183107 + 0.0im
    coeff3=-0.0 - 0.00674638241111651im
    memslots[4] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3]
    # Computing Bb4 = x*A+x*B2+x*B3
    coeff1=0.0 + 0.13340427306445612im
    coeff2=0.020226020298183107 + 0.0im
    coeff3=-0.0 - 0.00674638241111651im
    memslots[5] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3]
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[4],memslots[5])
    # Deallocating Ba4 in slot 4
    # Deallocating Bb4 in slot 5
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=2.6958430691533257 + 0.0im
    coeff2=0.0 + 0.05272871327381115im
    coeff3=-0.09896214548845832 + 0.0im
    coeff4=0.0 + 0.007295441446830946im
    coeff5=1.0 + 0.0im
    memslots[4] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[6]
    mul!(memslots[4],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=2.6958430691533257 + 0.0im
    coeff2=0.0 - 1.3591092616886926im
    coeff3=-0.09896214548845832 + 0.0im
    coeff4=0.0 + 0.015964794632994668im
    coeff5=1.0 + 0.0im
    # Smart lincomb recycle B4
    memslots[6] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[6]
    mul!(memslots[6],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[5],memslots[4],memslots[6])
    # Deallocating Ba5 in slot 4
    # Deallocating Bb5 in slot 6
    # Computing T2k7 = x*I+x*A+x*B2+x*B3+x*B5
    coeff1=-6.267569853502023 + 0.0im
    coeff2=0.0 + 2.521796947120981im
    coeff3=0.05786296656487002 + 0.0im
    coeff4=0.0 - 0.0776668640807187im
    coeff5=1.0 + 0.0im
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[5]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B5 in slot 5
    return memslots[1] # Returning T2k7
end

