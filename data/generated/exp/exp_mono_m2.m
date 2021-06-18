function output=exp_mono_m2(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: T2k2 B2 B3 T2k3 T2k5
    % Computing T2k2 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing B3 with operation: mult
    B3=B2 * A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.16666666666666666;
    T2k5= coeff1*T2k3+coeff2*B3;
    output=T2k5;
end

