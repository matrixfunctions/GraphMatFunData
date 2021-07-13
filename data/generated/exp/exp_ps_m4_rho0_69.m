function output=exp_ps_m4_rho0_69(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Bb4_2 Bb5_2 T2k2 B2 T2k3 Bb4_3 Bb5_3 B3 Bb4 B4 Bb5 B5 T2k7
    % Computing Bb4_2 with operation: lincomb
    coeff1=0.001388888888888889;
    coeff2=0.0001984126984126984;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.16666666666666666;
    coeff2=0.041666666666666664;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=2.48015873015873e-5;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.008333333333333333;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=B2 * A;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=2.7557319223985893e-6;
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

