function output=ps(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Bb4_2 Bb5_2 T2k2 B2 T2k3 Bb4_3 Bb5_3 B3 Bb4 B4 Bb5 B5 T2k7
    % Computing Bb4_2 with operation: lincomb
    coeff1=-0.0205078125;
    coeff2=0.01611328125;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.0625;
    coeff2=-0.0390625;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5;
    T2k2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.125;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.013092041015625;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02734375;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=B2 * A;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0109100341796875;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=B3 * Bb4;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Bb5= coeff1*Bb5_3+coeff2*B4;
    % Computing B5 with operation: mult
    B5=B3 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k7= coeff1*T2k3+coeff2*B5;
    output=T2k7;
end

