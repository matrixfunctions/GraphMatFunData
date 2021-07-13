function output=exp_sastre_m2(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: T2k2 Bb3_2 B2 Bb3 B3 T2k5
    % Computing T2k2 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k2= coeff1*I+coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.5;
    coeff2=0.16666666666666666;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.041666666666666664;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=B2 * Bb3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k5= coeff1*T2k2+coeff2*B3;
    output=T2k5;
end

