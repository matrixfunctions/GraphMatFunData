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

@inline function exp_sid_m5(A)
    return exp_sid_m5!(copy(A))
end

@inline function exp_sid_m5!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=8
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
    # Computation order: B2 B3 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 T2k8
    # Computing B2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing B3 with operation: mult
    mul!(memslots[3],memslots[1],memslots[2])
    # Computing Bb4 = x*A+x*B2+x*B3
    coeff1=5.374708803114821e-5
    coeff2=4.50085273957301e-6
    coeff3=1.16165883444488e-6
    memslots[4] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3]
    # Computing B4 with operation: mult
    mul!(memslots[5],memslots[3],memslots[4])
    # Deallocating Bb4 in slot 4
    # Computing Ba5 = x*A+x*B2+x*B3+x*B4
    coeff1=0.9418613214806352
    coeff2=0.06974348269544424
    coeff3=0.002005403977292901
    coeff4=1.0
    memslots[4] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3] .+ coeff4.*memslots[5]
    # Computing Bb5 = x*B2+x*B3+x*B4
    coeff1=-0.007544837153586671
    coeff2=0.002852960512714315
    coeff3=1.0
    memslots[6] .= coeff1.*memslots[2] .+ coeff2.*memslots[3] .+ coeff3.*memslots[5]
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[4],memslots[6])
    # Deallocating Ba5 in slot 4
    # Deallocating Bb5 in slot 6
    # Computing Ba6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.3112216227982407
    coeff2=0.0852839259083158
    coeff3=0.0292447258748138
    coeff4=1.829773504500424
    coeff5=1.0
    memslots[4] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3] .+ coeff4.*memslots[5] .+ coeff5.*memslots[7]
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.6865706355662834
    coeff2=0.1392249143769798
    coeff3=0.03151382711608315
    coeff4=11.173624766438472
    coeff5=1.0
    memslots[6] .= coeff1.*memslots[1] .+ coeff2.*memslots[2] .+ coeff3.*memslots[3] .+ coeff4.*memslots[5] .+ coeff5.*memslots[7]
    # Computing B6 with operation: mult
    mul!(memslots[8],memslots[4],memslots[6])
    # Deallocating Ba6 in slot 4
    # Deallocating Bb6 in slot 6
    # Computing T2k8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=1.0
    coeff2=1.0
    coeff3=0.2863243726334417
    coeff4=0.08776036732867759
    coeff5=0.18995526739487723
    coeff6=3.23337016308538
    coeff7=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[5] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 5
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    return memslots[1] # Returning T2k8
end

