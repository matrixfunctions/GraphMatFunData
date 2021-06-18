function output=bbcs_m1(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: T2k2 B2 T2k4
    % Computing T2k2 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.0 + 1i*-0.999999999976195;
    T2k2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.499999999992065 + 0.0im;
    T2k4= coeff1*T2k2+coeff2*B2;
    output=T2k4;
end

