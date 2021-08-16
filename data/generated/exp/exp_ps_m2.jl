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

@inline function exp_ps_m2(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_ps_m2!(A_copy)
end

@inline function exp_ps_m2!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
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
    # Computation order: A2 P1 C0 P0
    # Computing A2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing P1 = x*I+x*A+x*A2
    coeff1=0.5
    coeff2=0.16666666666666666
    coeff3=0.041666666666666664
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots2
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing C0 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating P1 in slot 3
    # Deallocating A2 in slot 2
    # Computing P0 = x*I+x*A+x*C0
    coeff1=1.0
    coeff2=1.0
    coeff3=1.0
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating C0 in slot 4
    return memslots1 # Returning P0
end

