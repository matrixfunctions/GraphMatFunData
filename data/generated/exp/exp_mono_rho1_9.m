function output=exp_mono_rho1_9(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: T2k2 B2 T2k3 B3 T2k4 B4 T2k5 B5 B6 T2k6 T2k8
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
    % Computing B3 with operation: mult
    B3=B2 * A;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.16666666666666666;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=B3 * A;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.041666666666666664;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=B4 * A;
    % Computing B6 with operation: mult
    B6=B5 * A;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.008333333333333333;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.001388888888888889;
    T2k8= coeff1*T2k6+coeff2*B6;
    output=T2k8;
end

