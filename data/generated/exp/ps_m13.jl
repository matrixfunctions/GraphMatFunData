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

@inline function ps_m13(A)
    return ps_m13!(copy(A))
end

@inline function ps_m13!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    max_memslots=10
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
    # Computation order: A2 A3 A4 A5 A6 A7 A8 P6 C5 P5 C4 P4 C3 P3 C2 P2 C1 P1 C0 P0
    # Computing A2 with operation: mult
    mul!(memslots[2],memslots[1],memslots[1])
    # Computing A3 with operation: mult
    mul!(memslots[3],memslots[1],memslots[2])
    # Computing A4 with operation: mult
    mul!(memslots[4],memslots[1],memslots[3])
    # Computing A5 with operation: mult
    mul!(memslots[5],memslots[1],memslots[4])
    # Computing A6 with operation: mult
    mul!(memslots[6],memslots[1],memslots[5])
    # Computing A7 with operation: mult
    mul!(memslots[7],memslots[1],memslots[6])
    # Computing A8 with operation: mult
    mul!(memslots[8],memslots[1],memslots[7])
    # Computing P6 = x*I+x*A+x*A2+x*A3+x*A4+x*A5+x*A6+x*A7+x*A8
    coeff1=8.055476070751236e-62
    coeff2=1.643974708316579e-63
    coeff3=3.287949416633158e-65
    coeff4=6.446959640457172e-67
    coeff5=1.2397999308571486e-68
    coeff6=2.3392451525606576e-70
    coeff7=4.331935467704922e-72
    coeff8=7.876246304918039e-74
    coeff9=1.4064725544496498e-75
    memslots[9] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4] .+ coeff6.*memslots[5] .+ coeff7.*memslots[6] .+ coeff8.*memslots[7] .+ coeff9.*memslots[8]
    mul!(memslots[9],true,I*coeff1,true,true)
    # Computing C5 with operation: mult
    mul!(memslots[10],memslots[9],memslots[8])
    # Deallocating P6 in slot 9
    # Computing P5 = x*I+x*A+x*A2+x*A3+x*A4+x*A5+x*A6+x*A7+x*C5
    coeff1=1.2256174391283858e-48
    coeff2=2.9893108271424046e-50
    coeff3=7.117406731291439e-52
    coeff4=1.6552108677421951e-53
    coeff5=3.7618428812322616e-55
    coeff6=8.359650847182804e-57
    coeff7=1.817315401561479e-58
    coeff8=3.866628513960594e-60
    coeff9=1.0
    # Smart lincomb recycle C5
    memslots[10] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4] .+ coeff6.*memslots[5] .+ coeff7.*memslots[6] .+ coeff8.*memslots[7] .+ coeff9.*memslots[10]
    mul!(memslots[10],true,I*coeff1,true,true)
    # Computing C4 with operation: mult
    mul!(memslots[9],memslots[10],memslots[8])
    # Deallocating P5 in slot 10
    # Computing P4 = x*I+x*A+x*A2+x*A3+x*A4+x*A5+x*A6+x*A7+x*C4
    coeff1=3.8003907548547434e-36
    coeff2=1.151633562077195e-37
    coeff3=3.387157535521162e-39
    coeff4=9.67759295863189e-41
    coeff5=2.6882202662866363e-42
    coeff6=7.265460179153071e-44
    coeff7=1.911963205040282e-45
    coeff8=4.902469756513544e-47
    coeff9=1.0
    # Smart lincomb recycle C4
    memslots[9] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4] .+ coeff6.*memslots[5] .+ coeff7.*memslots[6] .+ coeff8.*memslots[7] .+ coeff9.*memslots[9]
    mul!(memslots[9],true,I*coeff1,true,true)
    # Computing C3 with operation: mult
    mul!(memslots[10],memslots[9],memslots[8])
    # Deallocating P4 in slot 9
    # Computing P3 = x*I+x*A+x*A2+x*A3+x*A4+x*A5+x*A6+x*A7+x*C3
    coeff1=1.6117375710961184e-24
    coeff2=6.446950284384474e-26
    coeff3=2.4795962632247976e-27
    coeff4=9.183689863795546e-29
    coeff5=3.279889237069838e-30
    coeff6=1.1309962886447716e-31
    coeff7=3.7699876288159054e-33
    coeff8=1.216125041553518e-34
    coeff9=1.0
    # Smart lincomb recycle C3
    memslots[10] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4] .+ coeff6.*memslots[5] .+ coeff7.*memslots[6] .+ coeff8.*memslots[7] .+ coeff9.*memslots[10]
    mul!(memslots[10],true,I*coeff1,true,true)
    # Computing C2 with operation: mult
    mul!(memslots[9],memslots[10],memslots[8])
    # Deallocating P3 in slot 10
    # Computing P2 = x*I+x*A+x*A2+x*A3+x*A4+x*A5+x*A6+x*A7+x*C2
    coeff1=4.779477332387385e-14
    coeff2=2.8114572543455206e-15
    coeff3=1.5619206968586225e-16
    coeff4=8.22063524662433e-18
    coeff5=4.110317623312165e-19
    coeff6=1.9572941063391263e-20
    coeff7=8.896791392450574e-22
    coeff8=3.868170170630684e-23
    coeff9=1.0
    # Smart lincomb recycle C2
    memslots[9] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4] .+ coeff6.*memslots[5] .+ coeff7.*memslots[6] .+ coeff8.*memslots[7] .+ coeff9.*memslots[9]
    mul!(memslots[9],true,I*coeff1,true,true)
    # Computing C1 with operation: mult
    mul!(memslots[10],memslots[9],memslots[8])
    # Deallocating P2 in slot 9
    # Computing P1 = x*I+x*A+x*A2+x*A3+x*A4+x*A5+x*A6+x*A7+x*C1
    coeff1=2.48015873015873e-5
    coeff2=2.7557319223985893e-6
    coeff3=2.755731922398589e-7
    coeff4=2.505210838544172e-8
    coeff5=2.08767569878681e-9
    coeff6=1.6059043836821613e-10
    coeff7=1.1470745597729725e-11
    coeff8=7.647163731819816e-13
    coeff9=1.0
    # Smart lincomb recycle C1
    memslots[10] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4] .+ coeff6.*memslots[5] .+ coeff7.*memslots[6] .+ coeff8.*memslots[7] .+ coeff9.*memslots[10]
    mul!(memslots[10],true,I*coeff1,true,true)
    # Computing C0 with operation: mult
    mul!(memslots[9],memslots[10],memslots[8])
    # Deallocating P1 in slot 10
    # Deallocating A8 in slot 8
    # Computing P0 = x*I+x*A+x*A2+x*A3+x*A4+x*A5+x*A6+x*A7+x*C0
    coeff1=1.0
    coeff2=1.0
    coeff3=0.5
    coeff4=0.16666666666666666
    coeff5=0.041666666666666664
    coeff6=0.008333333333333333
    coeff7=0.001388888888888889
    coeff8=0.0001984126984126984
    coeff9=1.0
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[2] .+ coeff4.*memslots[3] .+ coeff5.*memslots[4] .+ coeff6.*memslots[5] .+ coeff7.*memslots[6] .+ coeff8.*memslots[7] .+ coeff9.*memslots[9]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating A2 in slot 2
    # Deallocating A3 in slot 3
    # Deallocating A4 in slot 4
    # Deallocating A5 in slot 5
    # Deallocating A6 in slot 6
    # Deallocating A7 in slot 7
    # Deallocating C0 in slot 9
    return memslots[1] # Returning P0
end

