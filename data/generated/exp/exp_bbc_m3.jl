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

@inline function exp_bbc_m3(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_bbc_m3!(A_copy)
end

@inline function exp_bbc_m3!(A)
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
    # Computation order: B2 Bb3 B3 Ba4 Bb4 B4 T2k6
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing Bb3 = x*A+x*B2
    coeff1=0.1083646567852278
    coeff2=0.02709116419630695
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing B3 with operation: mult
    mul!(memslots4,memslots2,memslots3)
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*B2+x*B3
    coeff1=0.6666666666666666
    coeff2=1.0
    memslots3 .= coeff1.*memslots2 .+ coeff2.*memslots4
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=0.546761457970724
    coeff2=0.16112557339541758
    coeff3=0.014090917158378206
    coeff4=0.0337927970108705
    # Smart lincomb recycle B3
    memslots4 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots4
    mul!(memslots4,true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots5,memslots3,memslots4)
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 4
    # Computing T2k6 = x*I+x*A+x*B2+x*B4
    coeff1=1.0
    coeff2=1.0
    coeff3=0.13549236135285073
    coeff4=1.0
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots5
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B4 in slot 5
    return memslots1 # Returning T2k6
end

