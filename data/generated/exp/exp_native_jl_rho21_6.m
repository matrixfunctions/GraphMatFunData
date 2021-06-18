function output=exp_native_jl_rho21_6(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: C A2 Ua2 Va2 A4 Va3 Ua3 Ub2 Vb2 A6 Ua Va Ub3 Ub Uc U Vb3 Vb V Z X P S1 S2
    % Computing C with operation: lincomb
    coeff1=0.25;
    coeff2=0.0;
    C= coeff1*A+coeff2*I;
    % Computing A2 with operation: mult
    A2=C * C;
    % Computing Ua2 with operation: lincomb
    coeff1=3.238237626624e16;
    coeff2=1.1873537964288e15;
    Ua2= coeff1*I+coeff2*A2;
    % Computing Va2 with operation: lincomb
    coeff1=6.476475253248e16;
    coeff2=7.7717703038976e15;
    Va2= coeff1*I+coeff2*A2;
    % Computing A4 with operation: mult
    A4=A2 * A2;
    % Computing Va3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.29060195264e14;
    Va3= coeff1*Va2+coeff2*A4;
    % Computing Ua3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.05594705216e13;
    Ua3= coeff1*Ua2+coeff2*A4;
    % Computing Ub2 with operation: lincomb
    coeff1=4.08408e7;
    coeff2=16380.0;
    Ub2= coeff1*A2+coeff2*A4;
    % Computing Vb2 with operation: lincomb
    coeff1=1.32324192e9;
    coeff2=960960.0;
    Vb2= coeff1*A2+coeff2*A4;
    % Computing A6 with operation: mult
    A6=A2 * A4;
    % Computing Ua with operation: lincomb
    coeff1=1.0;
    coeff2=3.352212864e10;
    Ua= coeff1*Ua3+coeff2*A6;
    % Computing Va with operation: lincomb
    coeff1=1.0;
    coeff2=6.704425728e11;
    Va= coeff1*Va3+coeff2*A6;
    % Computing Ub3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ub3= coeff1*Ub2+coeff2*A6;
    % Computing Ub with operation: mult
    Ub=Ub3 * A6;
    % Computing Uc with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Uc= coeff1*Ub+coeff2*Ua;
    % Computing U with operation: mult
    U=C * Uc;
    % Computing Vb3 with operation: lincomb
    coeff1=1.0;
    coeff2=182.0;
    Vb3= coeff1*Vb2+coeff2*A6;
    % Computing Vb with operation: mult
    Vb=Vb3 * A6;
    % Computing V with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    V= coeff1*Vb+coeff2*Va;
    % Computing Z with operation: lincomb
    coeff1=1.0;
    coeff2=-1.0;
    Z= coeff1*V+coeff2*U;
    % Computing X with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    X= coeff1*V+coeff2*U;
    % Computing P with operation: ldiv
    P=Z \ X;
    % Computing S1 with operation: mult
    S1=P * P;
    % Computing S2 with operation: mult
    S2=S1 * S1;
    output=S2;
end

