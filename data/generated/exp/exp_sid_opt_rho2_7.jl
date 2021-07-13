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

@inline function exp_sid_opt_rho2_7(A)
    return exp_sid_opt_rho2_7!(copy(A))
end

@inline function exp_sid_opt_rho2_7!(A)
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
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 Ba7 Bb7 B7 T2k9
    # Computing Bb2 = x*I+x*A
    coeff1=0.045435112444573475
    coeff2=1.0223246692300103
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=0.045435112444573475
    coeff2=1.0223246692300103
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=0.013779624734115105
    coeff2=1.0123617774427427
    coeff3=0.02096699267528694
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=0.0025872950992219526
    coeff2=0.013669635797719902
    coeff3=1.0123048626491173
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots[5],memslots[2],memslots[3])
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=0.0024179040551112
    coeff2=0.9995075315498182
    coeff3=0.0417319729524033
    coeff4=0.0010158820545039565
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=0.00023252325808871646
    coeff2=4.612679859250955e-5
    coeff3=0.010778337993824194
    coeff4=0.9999831553043256
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[2],memslots[3])
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=7.788069195751279e-6
    coeff2=8.668814613800501e-5
    coeff3=0.0005851748274979261
    coeff4=0.0003724118165732193
    coeff5=1.000002894879292
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=4.367862801332661e-5
    coeff2=4.771005078802816e-6
    coeff3=3.057699543578747e-7
    coeff4=1.4369102577037629e-8
    coeff5=2.8715019361588243e-10
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[2],memslots[3])
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-0.0014340288130826606
    coeff2=3.0214170368224798
    coeff3=0.28010545658135216
    coeff4=0.017678291021512645
    coeff5=0.0006690728448219196
    coeff6=0.7423612498678709
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.0009787233688326564
    coeff2=-0.096167679660017
    coeff3=-0.010643200832947241
    coeff4=-0.00024824310516588245
    coeff5=0.00012514140616432543
    coeff6=1.228607675289351
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[8],memslots[2],memslots[3])
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.0003158514159932758
    coeff2=8.290611574892452
    coeff3=1.7421134180734874
    coeff4=0.30262885192060646
    coeff5=0.02034868019711193
    coeff6=249.89631597753709
    coeff7=1.227713124647877
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-0.0026162370361472948
    coeff2=0.030992248965416677
    coeff3=0.009171582602004843
    coeff4=0.0014654665589331566
    coeff5=8.173851230520504e-5
    coeff6=1.119040499126371
    coeff7=1.4271746834360811e-5
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots[9],memslots[2],memslots[3])
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 3
    # Computing T2k9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.9996354364827498
    coeff2=1.0075667227584877
    coeff3=0.181554796011574
    coeff4=0.017961596596304457
    coeff5=0.0006512172198766862
    coeff6=6.761523475937711e-6
    coeff7=-0.036931956487199616
    coeff8=1.1186848822520254
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    # Deallocating B7 in slot 9
    return memslots[1] # Returning T2k9
end

