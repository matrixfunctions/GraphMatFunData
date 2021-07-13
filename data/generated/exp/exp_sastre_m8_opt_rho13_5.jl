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

@inline function exp_sastre_m8_opt_rho13_5(A)
    return exp_sastre_m8_opt_rho13_5!(copy(A))
end

@inline function exp_sastre_m8_opt_rho13_5!(A)
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
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 Ba7 Bb7 B7 Ba8 Bb8 B8 Ba9 Bb9 B9 T2k11
    # Computing Bb2 = x*I+x*A
    coeff1=-0.11490308438229017
    coeff2=0.15730697013710307
    memslots[2] .= coeff2.*memslots[1]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=-0.11490308438229017
    coeff2=0.15730697013710307
    memslots[3] .= coeff2.*memslots[1]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots[4],memslots[3],memslots[2])
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=0.004105729412589247
    coeff2=-0.006768382232040512
    coeff3=0.9999471251526917
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=-0.0028647512623116655
    coeff2=0.0009954035208072304
    coeff3=0.001211033078915847
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots[5],memslots[2],memslots[3])
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=0.027019231244596563
    coeff2=0.06170763799314927
    coeff3=-0.0018954361862412274
    coeff4=1.0067849094795742
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.013620620860785195
    coeff2=0.07638110624385147
    coeff3=0.14911649238384225
    coeff4=0.9051158545637227
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots[6],memslots[2],memslots[3])
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.0032451095696861268
    coeff2=-0.13204423967508136
    coeff3=0.9186197117187667
    coeff4=0.42211811124544796
    coeff5=-0.01259382817947474
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.011274427044379665
    coeff2=0.001517105603937254
    coeff3=0.001158435856997462
    coeff4=-0.0047467287622741695
    coeff5=-7.137183226929846e-5
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots[7],memslots[2],memslots[3])
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.011539085894536814
    coeff2=0.0935630221089867
    coeff3=0.06407916439324149
    coeff4=-0.285986613635287
    coeff5=0.044725605639523214
    coeff6=1.0008083736633122
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.02459649940876415
    coeff2=0.016824149148287185
    coeff3=0.03458106185410125
    coeff4=-0.047464943469750574
    coeff5=0.11471274073059659
    coeff6=1.0130277788963542
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots[8],memslots[2],memslots[3])
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=3.0826131069323353
    coeff2=0.2781411205676187
    coeff3=0.15593200742501928
    coeff4=7.925112972476668
    coeff5=1.001825018602058
    coeff6=-0.10421074361186333
    coeff7=0.16056529744528691
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.025970161938770917
    coeff2=0.05228210309031749
    coeff3=0.10608150003582509
    coeff4=0.049985371773829376
    coeff5=0.28766922644726883
    coeff6=1.936009578758007
    coeff7=1.0130307509989027
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots[9],memslots[2],memslots[3])
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 3
    # Computing Ba8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.9222548073365477
    coeff2=0.10238368416498037
    coeff3=0.1032554965257702
    coeff4=0.030058902112840666
    coeff5=-0.11384698074301966
    coeff6=0.011734571936531983
    coeff7=-0.05247830183077817
    coeff8=0.9972225209602774
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.9222548073365477
    coeff2=0.10238368416498037
    coeff3=0.1032554965257702
    coeff4=0.030058902112840666
    coeff5=-0.11384698074301966
    coeff6=0.011734571936531983
    coeff7=-0.05247830183077817
    coeff8=0.9972225209602774
    memslots[3] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9]
    mul!(memslots[3],true,I*coeff1,true,true)
    # Computing B8 with operation: mult
    mul!(memslots[10],memslots[2],memslots[3])
    # Deallocating Ba8 in slot 2
    # Deallocating Bb8 in slot 3
    # Computing Ba9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=0.3035081839334926
    coeff2=0.042128860424959606
    coeff3=0.005532874783353146
    coeff4=0.05011906670785025
    coeff5=0.03702751170872807
    coeff6=0.06620554678336489
    coeff7=0.15748603557212862
    coeff8=-0.03947847759964156
    coeff9=1.058382881818128
    memslots[2] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9] .+ coeff9.*memslots[10]
    mul!(memslots[2],true,I*coeff1,true,true)
    # Computing Bb9 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=-0.2738445106644693
    coeff2=-0.033002263370480726
    coeff3=0.004912635491924672
    coeff4=0.04449492500250622
    coeff5=0.03287724382394223
    coeff6=0.05879027030590609
    coeff7=0.13983408425497476
    coeff8=-0.03505336576885022
    coeff9=0.9397507089701836
    # Smart lincomb recycle B4
    memslots[6] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[6] .+ coeff6.*memslots[7] .+ coeff7.*memslots[8] .+ coeff8.*memslots[9] .+ coeff9.*memslots[10]
    mul!(memslots[6],true,I*coeff1,true,true)
    # Computing B9 with operation: mult
    mul!(memslots[3],memslots[2],memslots[6])
    # Deallocating Ba9 in slot 2
    # Deallocating Bb9 in slot 6
    # Computing T2k11 = x*I+x*A+x*B2+x*B3+x*B5+x*B6+x*B7+x*B8+x*B9
    coeff1=0.08072764899462631
    coeff2=0.02311006099515489
    coeff3=0.05527883142327172
    coeff4=-1.5378493201072977e-5
    coeff5=1.292606058383822e-5
    coeff6=6.578027159862427e-7
    coeff7=-1.3533497256574092e-7
    coeff8=1.9207658227240783e-6
    coeff9=0.9800057467529937
    # Smart lincomb recycle A
    memslots[1] .= coeff2.*memslots[1] .+ coeff3.*memslots[4] .+ coeff4.*memslots[5] .+ coeff5.*memslots[7] .+ coeff6.*memslots[8] .+ coeff7.*memslots[9] .+ coeff8.*memslots[10] .+ coeff9.*memslots[3]
    mul!(memslots[1],true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    # Deallocating B7 in slot 9
    # Deallocating B8 in slot 10
    # Deallocating B9 in slot 3
    return memslots[1] # Returning T2k11
end

