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

@inline function exp_sastre_m7_opt_rho6_0(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_sastre_m7_opt_rho6_0!(A_copy)
end

@inline function exp_sastre_m7_opt_rho6_0!(A)
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
    coeff1=0.0015638303136970531
    coeff2=0.5131399756353182
    memslots2 .= coeff2.*memslots1
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=0.0015638303136970531
    coeff2=0.5131399756353182
    memslots3 .= coeff2.*memslots1
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=0.0015678308682144713
    coeff2=-0.001644406122563135
    coeff3=1.0000648973491737
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=-0.009556418198624978
    coeff2=0.0008384641570218439
    coeff3=0.00012764063512709018
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots5,memslots2,memslots3)
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.0025650924538428376
    coeff2=0.0751072955748711
    coeff3=0.013533154589413957
    coeff4=1.1136335440836176
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=-0.029899374961962634
    coeff2=0.17661587215939123
    coeff3=0.05543256972765084
    coeff4=1.0216298885266226
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots6,memslots2,memslots3)
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.0021189804037908145
    coeff2=0.001037076938904245
    coeff3=1.0216723193069674
    coeff4=0.022323328476928236
    coeff5=-0.052946917852281214
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-0.016491731477218847
    coeff2=0.0008905549732892586
    coeff3=9.855583452895887e-5
    coeff4=-0.0008621015249925553
    coeff5=-2.1080367849877218e-5
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots7,memslots2,memslots3)
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.012857678730185047
    coeff2=0.20416569487077543
    coeff3=0.04205470017452324
    coeff4=0.01518879866885416
    coeff5=0.044060889132843686
    coeff6=1.0071150699762386
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.017396334722190973
    coeff2=0.04276963803433957
    coeff3=0.03219045127574837
    coeff4=-0.07421310609041179
    coeff5=0.09089674435798574
    coeff6=0.9929884164279147
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots8,memslots2,memslots3)
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing Ba7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=3.1087497397974926
    coeff2=0.3694214555352366
    coeff3=0.15935169296586307
    coeff4=7.9059413858030805
    coeff5=1.1026236327353052
    coeff6=-0.048941742425579186
    coeff7=0.24374173790191617
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb7 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.009133614148579657
    coeff2=0.14669027898511353
    coeff3=0.0775106388395176
    coeff4=0.020068924005133555
    coeff5=0.2544089375564568
    coeff6=1.9371031112856034
    coeff7=0.9819279402861139
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B7 with operation: mult
    mul!(memslots9,memslots2,memslots3)
    # Deallocating Ba7 in slot 2
    # Deallocating Bb7 in slot 3
    # Computing Ba8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.9717659119346435
    coeff2=0.03708072735488672
    coeff3=-0.0032737750163317585
    coeff4=0.021723450159739942
    coeff5=-0.13300038220488442
    coeff6=-0.005257171591711241
    coeff7=0.02345129553527665
    coeff8=0.9810542134903524
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=0.9717659119346435
    coeff2=0.03708072735488672
    coeff3=-0.0032737750163317585
    coeff4=0.021723450159739942
    coeff5=-0.13300038220488442
    coeff6=-0.005257171591711241
    coeff7=0.02345129553527665
    coeff8=0.9810542134903524
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8 .+ coeff8.*memslots9
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B8 with operation: mult
    mul!(memslots10,memslots2,memslots3)
    # Deallocating Ba8 in slot 2
    # Deallocating Bb8 in slot 3
    # Computing T2k10 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=-0.004180474678744686
    coeff2=-0.002570469257093357
    coeff3=-0.0009454274909569563
    coeff4=-0.003167848437529552
    coeff5=-0.003321896755863846
    coeff6=-0.015009539610819524
    coeff7=-0.02988218764945254
    coeff8=0.004141999495986352
    coeff9=1.002949767210698
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

