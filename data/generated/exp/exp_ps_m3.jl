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

@inline function exp_ps_m3(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_ps_m3!(A_copy)
end

@inline function exp_ps_m3!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=5
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: A2 A3 P1 C0 P0
    # Computing A2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing A3 with operation: mult
    mul!(memslots3,memslots1,memslots2)
    # Computing P1 = x*I+x*A+x*A2+x*A3
    coeff1=0.16666666666666666
    coeff2=0.041666666666666664
    coeff3=0.008333333333333333
    coeff4=0.001388888888888889
    memslots4 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3
    mul!(memslots4,true,I*coeff1,true,true)
    # Computing C0 with operation: mult
    mul!(memslots5,memslots4,memslots3)
    # Deallocating P1 in slot 4
    # Deallocating A3 in slot 3
    # Computing P0 = x*I+x*A+x*A2+x*C0
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=1.0
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots5
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating A2 in slot 2
    # Deallocating C0 in slot 5
    return memslots1 # Returning P0
end

