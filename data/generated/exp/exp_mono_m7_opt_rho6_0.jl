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

@inline function exp_mono_m7_opt_rho6_0(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_mono_m7_opt_rho6_0!(A_copy)
end

@inline function exp_mono_m7_opt_rho6_0!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=10
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    memslots7 = similar(A,T)
    memslots8 = similar(A,T)
    memslots9 = similar(A,T)
    memslots10 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 Ba7 Bb7 B7 Ba8 Bb8 B8 T2k10
    # Computing Bb2 = x*I+x*A
    coeff1=0.19853207576158147
    coeff2=0.4568156227748025
    memslots2 .= coeff2.*memslots1
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=0.19853207576158147
    coeff2=0.4568156227748025
    memslots3 .= coeff2.*memslots1
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=-0.03985250441225717
    coeff2=0.1604494657928954
    coeff3=0.9304976822330905
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=0.3044412670842123
    coeff2=0.4680830251883016
    coeff3=0.029166175784732254
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots5,memslots2,memslots3)
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.01206538482409284
    coeff2=-0.007524875465546422
    coeff3=0.32310560219905493
    coeff4=0.9639420128969182
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.011223283513386698
    coeff2=0.5004951105037614
    coeff3=0.0587301450366299
    coeff4=0.0020530677449521447
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots6,memslots2,memslots3)
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.0014350882470920096
    coeff2=-0.0029204003557249366
    coeff3=0.004322247579028231
    coeff4=0.013953282030764238
    coeff5=0.8220608344548307
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.8268548885673269
    coeff2=0.32692355515774696
    coeff3=0.05385438234818326
    coeff4=0.003348134917199928
    coeff5=4.728916919898747e-5
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots7,memslots2,memslots3)
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-0.00011974681141880265
    coeff2=-9.609608204465964e-5
    coeff3=0.004160420231960609
    coeff4=-0.13025928619739066
    coeff5=0.8471750409761895
    coeff6=0.8150363087806349
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-0.10687141220929966
    coeff2=0.49573353466641706
    coeff3=0.09669233528231447
    coeff4=0.008584243736221898
    coeff5=0.0002579652548561839
    coeff6=-6.616194274488254e-8
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots8,memslots2,memslots3)
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-4.539313763978297e-6
    coeff2=2.446796588580378e-6
    coeff3=4.215986081751763e-5
    coeff4=-0.0061270230013822065
    coeff5=0.06884543294512852
    coeff6=-0.1203307382956932
    coeff7=0.9948104309724103
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-8.127184169897754e-8
    coeff2=0.5009422778394859
    coeff3=0.10725741908100991
    coeff4=0.011203185515508135
    coeff5=0.0005125194848333935
    coeff6=1.4968901291255363e-5
    coeff7=1.941362084939232e-8
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots9,memslots2,memslots3)
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 3
    # Computing Ba8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.9834849370950302
    coeff2=0.40790543324377204
    coeff3=0.4202712650569836
    coeff4=0.15990304681757225
    coeff5=0.033732264254459236
    coeff6=0.008536371457893934
    coeff7=0.0010322498247622353
    coeff8=4.413245234287325e-5
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.9834849370950302
    coeff2=0.40790543324377204
    coeff3=0.4202712650569836
    coeff4=0.15990304681757225
    coeff5=0.033732264254459236
    coeff6=0.008536371457893934
    coeff7=0.0010322498247622353
    coeff8=4.413245234287325e-5
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B8 with operation: mult
    mul!(memslots10,memslots2,memslots3)
    # Deallocating Ba8 in slot 2
    # Deallocating Bb8 in slot 3
    # Computing T2k10 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=-0.00018687382272586462
    coeff2=-0.00010155826657765157
    coeff3=-0.0001809250904584009
    coeff4=-0.00011452178091543187
    coeff5=-3.544767352137871e-5
    coeff6=-1.1868405154805229e-5
    coeff7=-1.9506111159992416e-6
    coeff8=-1.0864858600469639e-7
    coeff9=1.0004223592531323
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9 .+ coeff9.*memslots10
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    # Deallocating B7 in slot 9
    # Deallocating B8 in slot 10
    return memslots1 # Returning T2k10
end

