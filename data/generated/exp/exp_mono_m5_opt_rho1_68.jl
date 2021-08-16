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

@inline function exp_mono_m5_opt_rho1_68(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_mono_m5_opt_rho1_68!(A_copy)
end

@inline function exp_mono_m5_opt_rho1_68!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=8
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    memslots7 = similar(A,T)
    memslots8 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    # Computation order: Bb2 Ba2 B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb6 B6 T2k8
    # Computing Bb2 = x*I+x*A
    coeff1=0.0010904030413732757
    coeff2=1.0008760081747117
    memslots2 .= coeff2.*memslots1
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Ba2 = x*I+x*A
    coeff1=0.0010904030413732757
    coeff2=1.0008760081747117
    memslots3 .= coeff2.*memslots1
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B2 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating Ba2 in slot 3
    # Deallocating Bb2 in slot 2
    # Computing Ba3 = x*I+x*A+x*B2
    coeff1=-4.6346414362462995e-5
    coeff2=0.0015183205977868746
    coeff3=0.9996379335010811
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb3 = x*I+x*A+x*B2
    coeff1=0.001530127478825033
    coeff2=0.9994108881899734
    coeff3=0.030492786462992225
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B3 with operation: mult
    mul!(memslots5,memslots2,memslots3)
    # Deallocating Ba3 in slot 2
    # Deallocating Bb3 in slot 3
    # Computing Ba4 = x*I+x*A+x*B2+x*B3
    coeff1=-5.439593588786742e-5
    coeff2=8.642508638357873e-5
    coeff3=0.0011897338053029579
    coeff4=0.9997653650686295
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb4 = x*I+x*A+x*B2+x*B3
    coeff1=0.0010337804923906489
    coeff2=0.9970317188215052
    coeff3=0.08229676673753394
    coeff4=0.0031469743629754005
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B4 with operation: mult
    mul!(memslots6,memslots2,memslots3)
    # Deallocating Ba4 in slot 2
    # Deallocating Bb4 in slot 3
    # Computing Ba5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=-8.30025660582309e-6
    coeff2=-1.614854194638885e-5
    coeff3=-1.2536139794669188e-5
    coeff4=0.0010411131935665154
    coeff5=0.9999723778369385
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb5 = x*I+x*A+x*B2+x*B3+x*B4
    coeff1=0.0007031814361206821
    coeff2=1.0010882362282207
    coeff3=0.05182299699949662
    coeff4=0.007748443487686225
    coeff5=0.0002736987094186782
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B5 with operation: mult
    mul!(memslots7,memslots2,memslots3)
    # Deallocating Ba5 in slot 2
    # Deallocating Bb5 in slot 3
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-1.2408057646229667e-6
    coeff2=-3.1395014660350098e-6
    coeff3=-1.8699534803824349e-6
    coeff4=2.3223775362920453e-6
    coeff5=0.0007434150049329373
    coeff6=0.999999459268548
    memslots2 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots2,true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-2.7188726824232142e-6
    coeff2=0.9998021510172249
    coeff3=0.06065097803531471
    coeff4=0.005091212395247093
    coeff5=0.00034385565143820764
    coeff6=2.2549509122548797e-6
    memslots3 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7
    mul!(memslots3,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots8,memslots2,memslots3)
    # Deallocating Ba6 in slot 2
    # Deallocating Bb6 in slot 3
    # Computing T2k8 = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=0.9999994211561953
    coeff2=0.9989210320213427
    coeff3=0.49825582022614506
    coeff4=0.16623132434744778
    coeff5=0.036597641936229344
    coeff6=0.004190644602100282
    coeff7=0.0004888481055895867
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots4 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots8
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 4
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 8
    return memslots1 # Returning T2k8
end

