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

@inline function sastre_m3(A)
    return sastre_m3!(copy(A))
end

@inline function sastre_m3!(A)
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
    # Computation order: B2 Bb3 B3 Ba4 Bb4 B4 T2k6
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing Bb3 = x*A+x*B2
    coeff1=0.01992047682223989
    coeff2=0.004980119205559973
    memslots[3] .= coeff1.*memslots[1] .+ coeff2.*memslots[2]
    # Computing B3 with operation: mult
    mul!(memslots[4],memslots[2],memslots[3])
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*A+x*B2+x*B3
    coeff1=0.8765009801785554
    coeff2=0.07665265321119147
    coeff3=1.0
    memslots[3] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[4]
    # Computing Bb4 = x*B2+x*B3
    coeff1=0.1225521150112075
    coeff2=1.0
    memslots[5] .= coeff1.*memslots[2] .+ coeff2.*memslots[4]
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[3],memslots[5])
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 5
    # Computing T2k6 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=2.974307204847627
    coeff5=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[4] .+ coeff5.*memslots[6]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 4
    # Deallocating B4 in slot 6
    return memslots[1] # Returning T2k6
end

