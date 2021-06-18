function output=ps_m4(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: B_0_1 B_1_1 B_2_1 A2 B_0_2 B_2_2 B_1_2 A3 P2 C1 P1 C0 P0
    % Computing B_0_1 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    B_0_1= coeff1*I+coeff2*A;
    % Computing B_1_1 with operation: lincomb
    coeff1=0.16666666666666666;
    coeff2=0.041666666666666664;
    B_1_1= coeff1*I+coeff2*A;
    % Computing B_2_1 with operation: lincomb
    coeff1=0.001388888888888889;
    coeff2=0.0001984126984126984;
    B_2_1= coeff1*I+coeff2*A;
    % Computing A2 with operation: mult
    A2=A * A;
    % Computing B_0_2 with operation: lincomb
    coeff1=0.5;
    coeff2=1.0;
    B_0_2= coeff1*A2+coeff2*B_0_1;
    % Computing B_2_2 with operation: lincomb
    coeff1=2.48015873015873e-5;
    coeff2=1.0;
    B_2_2= coeff1*A2+coeff2*B_2_1;
    % Computing B_1_2 with operation: lincomb
    coeff1=0.008333333333333333;
    coeff2=1.0;
    B_1_2= coeff1*A2+coeff2*B_1_1;
    % Computing A3 with operation: mult
    A3=A * A2;
    % Computing P2 with operation: lincomb
    coeff1=2.7557319223985893e-6;
    coeff2=1.0;
    P2= coeff1*A3+coeff2*B_2_2;
    % Computing C1 with operation: mult
    C1=P2 * A3;
    % Computing P1 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    P1= coeff1*C1+coeff2*B_1_2;
    % Computing C0 with operation: mult
    C0=P1 * A3;
    % Computing P0 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    P0= coeff1*C0+coeff2*B_0_2;
    output=P0;
end

