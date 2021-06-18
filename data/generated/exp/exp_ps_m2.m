function output=exp_ps_m2(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: B_0_1 B_1_1 A2 P1 C0 P0
    % Computing B_0_1 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    B_0_1= coeff1*I+coeff2*A;
    % Computing B_1_1 with operation: lincomb
    coeff1=0.5;
    coeff2=0.16666666666666666;
    B_1_1= coeff1*I+coeff2*A;
    % Computing A2 with operation: mult
    A2=A * A;
    % Computing P1 with operation: lincomb
    coeff1=0.041666666666666664;
    coeff2=1.0;
    P1= coeff1*A2+coeff2*B_1_1;
    % Computing C0 with operation: mult
    C0=P1 * A2;
    % Computing P0 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    P0= coeff1*C0+coeff2*B_0_1;
    output=P0;
end

