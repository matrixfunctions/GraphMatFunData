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

@inline function denman_beavers(A)
    return denman_beavers!(copy(A))
end

@inline function denman_beavers!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=4
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
    # Computation order: A_shift Xinv0 Y1 X1 Yinv1 Xinv1 Y2 Yinv2 X3
    # Computing A_shift = x*A+x*I
    coeff1=1.0
    coeff2=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff1.*memslots[1]
    mul!(memslots[1],true,I*coeff2,true,true)
    # Computing Xinv0 with operation: ldiv
    memslots[2]=inv(memslots[1])
    # Computing Y1 = x*I+x*Xinv0
    coeff1=0.5
    coeff2=0.5
    # Smart lincomb recycle Xinv0
    memslots[2] .= coeff2.*memslots[2]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing X1 = x*A_shift+x*I
    coeff1=0.5
    coeff2=0.5
    # Smart lincomb recycle A_shift
    memslots[1] .= coeff1.*memslots[1]
    mul!(memslots[1],true,I*coeff2,true,true)
    # Computing Yinv1 with operation: ldiv
    memslots[3]=inv(memslots[2])
    # Computing Xinv1 with operation: ldiv
    memslots[4]=inv(memslots[1])
    # Computing Y2 = x*Y1+x*Xinv1
    coeff1=0.5
    coeff2=0.5
    # Smart lincomb recycle Y1
    memslots[2] .= coeff1.*memslots[2] .+ coeff2.*memslots[4]
    # Deallocating Xinv1 in slot 4
    # Computing Yinv2 with operation: ldiv
    memslots[4]=inv(memslots[2])
    # Deallocating Y2 in slot 2
    # Computing X3 = x*X1+x*Yinv1+x*Yinv2
    coeff1=0.25
    coeff2=0.25
    coeff3=0.5
    # Smart lincomb recycle X1
    memslots[1] .= coeff1.*memslots[1] .+ coeff2.*memslots[3] .+ coeff3.*memslots[4]
    # Deallocating Yinv1 in slot 3
    # Deallocating Yinv2 in slot 4
    return memslots[1] # Returning X3
end

