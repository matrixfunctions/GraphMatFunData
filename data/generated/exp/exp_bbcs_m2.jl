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

@inline function exp_bbcs_m2(A)
    T=promote_type(eltype(A),ComplexF64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_bbcs_m2!(A_copy)
end

@inline function exp_bbcs_m2!(A)
    T=promote_type(eltype(A),ComplexF64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=4
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: B2 Bb3 B3 T2k5
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing Bb3 = x*A+x*B2
    coeff1=0.0 + 0.16666657785001893im
    coeff2=0.04166664890333649 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing B3 with operation: mult
    mul!(memslots4,memslots2,memslots3)
    # Deallocating Bb3 in slot 3
    # Computing T2k5 = x*I+x*A+x*B2+x*B3
    coeff1=1.0 + 0.0im
    coeff2=0.0 - 0.9999999999998107im
    coeff3=-0.4999999999999432 + 0.0im
    coeff4=1.0 + 0.0im
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots4
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 4
    return memslots1 # Returning T2k5
end

