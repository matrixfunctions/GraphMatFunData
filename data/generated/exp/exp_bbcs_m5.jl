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

@inline function exp_bbcs_m5(A)
    T=promote_type(eltype(A),ComplexF64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_bbcs_m5!(A_copy)
end

@inline function exp_bbcs_m5!(A)
    T=promote_type(eltype(A),ComplexF64) # Make it work for many 'bigger' types (matrices and scalars)
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
    # Computation order: B2 B3 Ba5_4 B4 Bb5 B5 Ba6 Bb6 B6 T2k8
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing B3 with operation: mult
    mul!(memslots3,memslots2,memslots1)
    # Computing Ba5_4 = x*A+x*B2+x*B3
    coeff1=0.12 + 0.0im
    coeff2=0.0 - 0.008774760968797039im
    coeff3=-0.0009784845352378095 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3
    # Computing B4 with operation: mult
    mul!(memslots5,memslots3,memslots3)
    # Computing Bb5 = x*B2+x*B3+x*B4
    coeff1=0.0 - 0.12395369585828313im
    coeff2=-0.011202694841085593 + 0.0im
    coeff3=-0.0 - 1.2367240538259895e-5im
    memslots6 .= coeff1.*memslots2 .+ coeff2.*memslots3 .+ coeff3.*memslots5
    # Computing B5 with operation: mult
    mul!(memslots7,memslots4,memslots6)
    # Deallocating Ba5_4 in slot 4
    # Deallocating Bb5 in slot 6
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.3420232802536553 + 0.0im
    coeff2=0.0 - 0.2851997796332415im
    coeff3=0.047347067331271094 + 0.0im
    coeff4=0.0 - 0.022186600635366212im
    coeff5=-9.74758985615379e-6 + 0.0im
    coeff6=1.0 + 0.0im
    memslots4 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots5 .+ coeff6.*memslots7
    mul!(memslots4,true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=2.9237775839655367 + 0.0im
    coeff2=0.0 + 1.4451330034748826im
    coeff3=0.1240818356655045 + 0.0im
    coeff4=0.0 - 0.019571570936427238im
    coeff5=2.425253007433925e-5 + 0.0im
    coeff6=1.0 + 0.0im
    # Smart lincomb recycle B5
    memslots7 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots5 .+ coeff6.*memslots7
    mul!(memslots7,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots6,memslots4,memslots7)
    # Deallocating Ba6 in slot 4
    # Deallocating Bb6 in slot 7
    # Computing T2k8 = x*A+x*B2+x*B3+x*B4+x*B6
    coeff1=0.0 - 0.6604084076077131im
    coeff2=-1.093022784715649 + 0.0im
    coeff3=0.0 + 0.2537715581771087im
    coeff4=0.0005437426743473122 + 0.0im
    coeff5=1.0 + 0.0im
    # Smart lincomb recycle A
    memslots1 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 5
    # Deallocating B6 in slot 6
    return memslots1 # Returning T2k8
end

