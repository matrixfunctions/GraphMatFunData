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

@inline function exp_bbc_m5(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); A_copy .= A;
    return exp_bbc_m5!(A_copy)
end

@inline function exp_bbc_m5!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=7
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    memslots7 = similar(A,T)
    # Assign precomputed nodes memslots 
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    value_one=ValueOne()
    # Computation order: B2 B3 Ba5_4 B4 Bb5 B5 Ba6 Bb6 B6 T2k8
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing B3 with operation: mult
    mul!(memslots3,memslots1,memslots2)
    # Computing Ba5_4 = x*A+x*B2+x*B3
    coeff1=-0.10036558103014462
    coeff2=-0.00802924648241157
    coeff3=-0.00089213849804573
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3
    # Computing B4 with operation: mult
    mul!(memslots5,memslots3,memslots3)
    # Computing Bb5 = x*B2+x*B3+x*B4
    coeff1=-0.09233646193671186
    coeff2=-0.016936493900208172
    coeff3=-1.400867981820361e-5
    memslots6 .= coeff1.*memslots2 .+ coeff2.*memslots3 .+ coeff3.*memslots5
    # Computing B5 with operation: mult
    mul!(memslots7,memslots4,memslots6)
    # Deallocating Ba5_4 in slot 4
    # Deallocating Bb5 in slot 6
    # Computing Ba6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-11.058071288535288
    coeff2=1.6125176868819238
    coeff3=0.12477411482493252
    coeff4=0.02257315581805103
    coeff5=1.9579475957000978e-5
    coeff6=1.0
    memslots4 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots5 .+ coeff6.*memslots7
    mul!(memslots4,true,I*coeff1,true,true)
    # Computing Bb6 = x*I+x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-0.09043168323908106
    coeff2=-0.06764045190713819
    coeff3=0.06759613017704597
    coeff4=0.029555257042931552
    coeff5=-1.391802575160607e-5
    coeff6=1.0
    # Smart lincomb recycle B5
    memslots7 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots3 .+ coeff5.*memslots5 .+ coeff6.*memslots7
    mul!(memslots7,true,I*coeff1,true,true)
    # Computing B6 with operation: mult
    mul!(memslots6,memslots4,memslots7)
    # Deallocating Ba6 in slot 4
    # Deallocating Bb6 in slot 7
    # Computing T2k8 = x*A+x*B2+x*B3+x*B4+x*B6
    coeff1=0.3978497494996451
    coeff2=1.3678377846041172
    coeff3=0.49828962252538267
    coeff4=-0.0006378981945947233
    coeff5=1.0
    # Smart lincomb recycle A
    memslots1 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots5 .+ coeff5.*memslots6
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 3
    # Deallocating B4 in slot 5
    # Deallocating B6 in slot 6
    return memslots1 # Returning T2k8
end

