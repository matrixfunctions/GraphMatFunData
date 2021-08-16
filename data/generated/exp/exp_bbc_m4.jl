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

@inline function exp_bbc_m4(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_bbc_m4!(A_copy)
end

@inline function exp_bbc_m4!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=6
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: B2 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 T2k7
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing B3 with operation: mult
    mul!(memslots3,memslots2,memslots1)
    # Computing Ba4 = x*A+x*B2+x*B3
    coeff1=-0.13181061013830184
    coeff2=-0.02027855540589259
    coeff3=-0.006759518468630863
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3
    # Computing Bb4 = x*A+x*B2+x*B3
    coeff1=-0.13181061013830184
    coeff2=-0.02027855540589259
    coeff3=-0.006759518468630863
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3
    # Computing B4 with operation: mult
    mul!(memslots6,memslots4,memslots5)
    # Deallocating Ba4 in slot 4
    # Deallocating Bb4 in slot 5
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=4.811693118299809
    coeff2=1.1510994882542136
    coeff3=0.03318960838392776
    coeff4=0.012516177931579242
    coeff5=1.0
    memslots4 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots6
    mul!(memslots4,true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.21169311829980944
    coeff2=0.15822438471572672
    coeff3=0.1656351694367274
    coeff4=0.010786277931579243
    coeff5=1.0
    # Smart lincomb recycle B4
    memslots6 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots6
    mul!(memslots6,true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots5,memslots4,memslots6)
    # Deallocating Ba5 in slot 4
    # Deallocating Bb5 in slot 6
    # Computing T2k7 = x*I+x*A+x*B2+x*B3+x*B5
    coeff1=-0.018602320514620553
    coeff2=-0.005007023225733177
    coeff3=-0.5734201229605222
    coeff4=-0.13339969394389206
    coeff5=1.0
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots5
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B5 in slot 5
    return memslots1 # Returning T2k7
end

