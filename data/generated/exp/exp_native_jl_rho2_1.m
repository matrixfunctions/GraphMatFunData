function output=exp_native_jl_rho2_1(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: A2 Ua2 V2 A4 Ua3 V3 A6 Ua4 V4 A8 Ua U V Z X P
    % Computing A2 with operation: mult
    A2=A * A;
    % Computing Ua2 with operation: lincomb
    coeff1=8.8216128e9;
    coeff2=3.027024e8;
    Ua2= coeff1*I+coeff2*A2;
    % Computing V2 with operation: lincomb
    coeff1=1.76432256e10;
    coeff2=2.0756736e9;
    V2= coeff1*I+coeff2*A2;
    % Computing A4 with operation: mult
    A4=A2 * A2;
    % Computing Ua3 with operation: lincomb
    coeff1=1.0;
    coeff2=2.16216e6;
    Ua3= coeff1*Ua2+coeff2*A4;
    % Computing V3 with operation: lincomb
    coeff1=1.0;
    coeff2=3.027024e7;
    V3= coeff1*V2+coeff2*A4;
    % Computing A6 with operation: mult
    A6=A2 * A4;
    % Computing Ua4 with operation: lincomb
    coeff1=1.0;
    coeff2=3960.0;
    Ua4= coeff1*Ua3+coeff2*A6;
    % Computing V4 with operation: lincomb
    coeff1=1.0;
    coeff2=110880.0;
    V4= coeff1*V3+coeff2*A6;
    % Computing A8 with operation: mult
    A8=A2 * A6;
    % Computing Ua with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ua= coeff1*Ua4+coeff2*A8;
    % Computing U with operation: mult
    U=Ua * A;
    % Computing V with operation: lincomb
    coeff1=1.0;
    coeff2=90.0;
    V= coeff1*V4+coeff2*A8;
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
    output=P;
end

