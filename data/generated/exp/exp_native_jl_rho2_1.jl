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

@inline function exp_native_jl_rho2_1(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_native_jl_rho2_1!(A_copy)
end

@inline function exp_native_jl_rho2_1!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=7
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    memslots7 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: A2 A4 A6 A8 Ua U V Z X P
    # Computing A2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing A4 with operation: mult
    mul!(memslots3,memslots2,memslots2)
    # Computing A6 with operation: mult
    mul!(memslots4,memslots2,memslots3)
    # Computing A8 with operation: mult
    mul!(memslots5,memslots2,memslots4)
    # Computing Ua = x*I+x*A2+x*A4+x*A6+x*A8
    coeff1=8.8216128e9
    coeff2=3.027024e8
    coeff3=2.16216e6
    coeff4=3960.0
    coeff5=1.0
    memslots6 .= coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4 .+ coeff5.*memslots5
    mul!(memslots6,true,I*coeff1,true,true)
    # Computing U with operation: mult
    mul!(memslots7,memslots6,memslots1)
    # Deallocating Ua in slot 6
    # Deallocating A in slot 1
    # Computing V = x*I+x*A2+x*A4+x*A6+x*A8
    coeff1=1.76432256e10
    coeff2=2.0756736e9
    coeff3=3.027024e7
    coeff4=110880.0
    coeff5=90.0
    # Smart lincomb recycle A2
    memslots2 .= coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4 .+ coeff5.*memslots5
    mul!(memslots2,true,I*coeff1,true,true)
    # Deallocating A4 in slot 3
    # Deallocating A6 in slot 4
    # Deallocating A8 in slot 5
    # Computing Z = x*V+x*U
    coeff1=1.0
    coeff2=-1.0
    memslots1 .= coeff1.*memslots2 .+ coeff2.*memslots7
    # Computing X = x*V+x*U
    coeff1=1.0
    coeff2=1.0
    # Smart lincomb recycle V
    memslots2 .= coeff1.*memslots2 .+ coeff2.*memslots7
    # Deallocating U in slot 7
    # Computing P with operation: ldiv
    memslots3 .=memslots1\memslots2
    # Deallocating Z in slot 1
    # Deallocating X in slot 2
    return memslots3 # Returning P
end

