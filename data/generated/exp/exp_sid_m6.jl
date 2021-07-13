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

@inline function exp_sid_m6(A)
    return exp_sid_m6!(copy(A))
end

@inline function exp_sid_m6!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=9
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
    # Computation order: B2 B3 B4 Bb5 B5 Ba6 Bb7_6 Bb6 B6 Ba7 B7 T2k9
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing B3 with operation: mult
    mul!(memslots[3],memslots[1],memslots[2])
    # Computing B4 with operation: mult
    mul!(memslots[4],memslots[1],memslots[3])
    # Computing Bb5 = x*A+x*B2+x*B3+x*B4
    coeff1=2.294895435403922e-5
    coeff2=1.406952242413849e-6
    coeff3=9.379681616092325e-8
    coeff4=1.172460202011541e-8
    memslots[5] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3] .+ coeff4.*memslots[4]
    # Computing B5 with operation: mult
    mul!(memslots[6],memslots[4],memslots[5])
    # Deallocating Bb5 in slot 5
    # Computing Ba6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=2.865001388641538
    coeff2=0.1952545843107103
    coeff3=0.01430688980356062
    coeff4=0.002024281516007681
    coeff5=1.0
    memslots[5] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3] .+ coeff4.*memslots[4] .+ coeff5.*memslots[6]
    # Computing Bb7_6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.03655234395347475
    coeff2=0.01606091400855144
    coeff3=0.0001758035313846159
    coeff4=0.0002919349464582001
    coeff5=1.0
    memslots[7] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3] .+ coeff4.*memslots[4] .+ coeff5.*memslots[6]
    # Computing Bb6 = x*B2+x*B3+x*B4+x*B5
    coeff1=0.02721930992200371
    coeff2=0.002547056607231984
    coeff3=-0.001204349003694297
    coeff4=1.0
    memslots[8] .= coeff1.*memslots[2] .+ coeff2.*memslots[3] .+ coeff3.*memslots[4] .+ coeff4.*memslots[6]
    # Computing B6 with operation: mult
    mul!(memslots[9],memslots[5],memslots[8])
    # Deallocating Ba6 in slot 5
    # Deallocating Bb6 in slot 8
    # Computing Ba7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=8.29008575139441
    coeff2=1.739158441630994
    coeff3=0.1965098904519709
    coeff4=0.02018492049443954
    coeff5=249.896909254999
    coeff6=1.0
    # Smart lincomb recycle B5
    memslots[6] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3] .+ coeff4.*memslots[4] .+ coeff5.*memslots[6] .+ coeff6.*memslots[9]
    # Deallocating B6 in slot 9
    # Computing B7 with operation: mult
    mul!(memslots[5],memslots[6],memslots[7])
    # Deallocating Ba7 in slot 6
    # Deallocating Bb7_6 in slot 7
    # Computing T2k9 = x*I+x*A+x*B2+x*B3+x*B4+x*B7
    coeff1=1.0
    coeff2=1.0
    coeff3=0.1969779342112314
    coeff4=-0.03005000525808178
    coeff5=0.002243394407902074
    coeff6=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4] .+ coeff6.*memslots[5]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 4
    # Deallocating B7 in slot 5
    return memslots[1] # Returning T2k9
end

