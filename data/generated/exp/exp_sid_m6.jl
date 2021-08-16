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
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_sid_m6!(A_copy)
end

@inline function exp_sid_m6!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=9
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
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: B2 B3 B4 Bb5 B5 Ba6 Bb7_6 Bb6 B6 Ba7 B7 T2k9
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing B3 with operation: mult
    mul!(memslots3,memslots1,memslots2)
    # Computing B4 with operation: mult
    mul!(memslots4,memslots1,memslots3)
    # Computing Bb5 = x*A+x*B2+x*B3+x*B4
    coeff1=2.294895435403922e-5
    coeff2=1.406952242413849e-6
    coeff3=9.379681616092325e-8
    coeff4=1.172460202011541e-8
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4
    # Computing B5 with operation: mult
    mul!(memslots6,memslots4,memslots5)
    # Deallocating Bb5 in slot 5
    # Computing Ba6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=2.865001388641538
    coeff2=0.1952545843107103
    coeff3=0.01430688980356062
    coeff4=0.002024281516007681
    coeff5=1.0
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4 .+ coeff5.*memslots6
    # Computing Bb7_6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=0.03655234395347475
    coeff2=0.01606091400855144
    coeff3=0.0001758035313846159
    coeff4=0.0002919349464582001
    coeff5=1.0
    memslots7 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4 .+ coeff5.*memslots6
    # Computing Bb6 = x*B2+x*B3+x*B4+x*B5
    coeff1=0.02721930992200371
    coeff2=0.002547056607231984
    coeff3=-0.001204349003694297
    coeff4=1.0
    memslots8 .= coeff1.*memslots2 .+ coeff2.*memslots3 .+ coeff3.*memslots4 .+ coeff4.*memslots6
    # Computing B6 with operation: mult
    mul!(memslots9,memslots5,memslots8)
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
    memslots6 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4 .+ coeff5.*memslots6 .+ coeff6.*memslots9
    # Deallocating B6 in slot 9
    # Computing B7 with operation: mult
    mul!(memslots5,memslots6,memslots7)
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
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots4 .+ coeff6.*memslots5
    mul!(memslots1,true,I*coeff1,true,true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 4
    # Deallocating B7 in slot 5
    return memslots1 # Returning T2k9
end

