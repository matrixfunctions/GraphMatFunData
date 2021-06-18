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

@inline function bbcs_m7(A)
    return bbcs_m7!(copy(A))
end

@inline function bbcs_m7!(A)
    T=promote_type(eltype(A),ComplexF64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=7
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
    # Computation order: Bb2 Ba2 B2 Bb3 B3 Ba5_4 B4 Bb5 B5 Ba6 Bb6 B6 Ba7 Bb7 B7 B8 T2k10
    # Computing Bb2 = x*I+x*A
    coeff1=0.0 + 0.0im
    coeff2=0.25 + 0.0im
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=0.0 + 0.0im
    coeff2=0.25 + 0.0im
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Bb3 = x*A+x*B2
    coeff1=0.25 + 0.0im
    coeff2=0.0 + 0.0im
    memslots[2] .= coeff1.*memslots[1] .+ coeff2.*memslots[4]
    # Computing B3 with operation: mult
    mul!(memslots[3],memslots[4],memslots[2])
    # Deallocating Bb3 in slot 2
    # Computing Ba5_4 = x*A+x*B2+x*B3
    coeff1=0.03 + 0.0im
    coeff2=0.0 - 0.008774760968797039im
    coeff3=-0.0009784845352378095 + 0.0im
    memslots[2] .= coeff1.*memslots[1] .+ coeff2.*memslots[4] .+ coeff3.*memslots[3]
    # Computing B4 with operation: mult
    mul!(memslots[5],memslots[3],memslots[3])
    # Computing Bb5 = x*B2+x*B3+x*B4
    coeff1=0.0 - 0.12395369585828313im
    coeff2=-0.011202694841085593 + 0.0im
    coeff3=-0.0 - 1.2367240538259895e-5im
    memslots[6] .= coeff1.*memslots[4] .+ coeff2.*memslots[3] .+ coeff3.*memslots[5]
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[2],memslots[6])
    # Deallocating Ba5_4 in slot 2
    # Deallocating Bb5 in slot 6
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.3420232802536553 + 0.0im
    coeff2=0.0 - 0.07129994490831038im
    coeff3=0.047347067331271094 + 0.0im
    coeff4=0.0 - 0.022186600635366212im
    coeff5=-9.74758985615379e-6 + 0.0im
    coeff6=1.0 + 0.0im
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[3] .+ coeff5.*memslots[5] .+ coeff6.*memslots[7]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=2.9237775839655367 + 0.0im
    coeff2=0.0 + 0.36128325086872065im
    coeff3=0.1240818356655045 + 0.0im
    coeff4=0.0 - 0.019571570936427238im
    coeff5=2.425253007433925e-5 + 0.0im
    coeff6=1.0 + 0.0im
    # Smart lincomb recycle B5
    memslots[7] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[3] .+ coeff5.*memslots[5] .+ coeff6.*memslots[7]
    mul!(memslots[7],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[6],memslots[2],memslots[7])
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 7
    # Computing Ba7 = x*A+x*B2+x*B3+x*B4+x*B6
    coeff1=0.0 - 0.16510210190192828im
    coeff2=-1.093022784715649 + 0.0im
    coeff3=0.0 + 0.2537715581771087im
    coeff4=0.0005437426743473122 + 0.0im
    coeff5=1.0 + 0.0im
    memslots[2] .= coeff1.*memslots[1] .+ coeff2.*memslots[4] .+ coeff3.*memslots[3] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    # Computing Bb7 = x*A+x*B2+x*B3+x*B4+x*B6
    coeff1=0.0 - 0.16510210190192828im
    coeff2=-1.093022784715649 + 0.0im
    coeff3=0.0 + 0.2537715581771087im
    coeff4=0.0005437426743473122 + 0.0im
    coeff5=1.0 + 0.0im
    # Smart lincomb recycle A
    memslots[1] .= coeff1.*memslots[1] .+ coeff2.*memslots[4] .+ coeff3.*memslots[3] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 5
    # Deallocating B6 in slot 6
    # Computing B7 with operation: mult
    mul!(memslots[3],memslots[2],memslots[1])
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 1
    # Computing B8 with operation: mult
    mul!(memslots[1],memslots[3],memslots[3])
    # Computing T2k10 = x*B7+x*B8
    coeff1=0.0 + 0.0im
    coeff2=1.0 + 0.0im
    # Smart lincomb recycle B7
    memslots[3] .= coeff1.*memslots[3] .+ coeff2.*memslots[1]
    # Deallocating B8 in slot 1
    return memslots[3] # Returning T2k10
end

