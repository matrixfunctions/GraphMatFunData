function output=exp_bbcs_m2(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: T2k2 B2 Bb3 B3 T2k3 T2k5
    % Computing T2k2 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.0 + 1i*-0.9999999999998107;
    T2k2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing Bb3 with operation: lincomb
    coeff1=0.0 + 1i*0.16666657785001893;
    coeff2=0.04166664890333649 + 0.0im;
    Bb3= coeff1*A+coeff2*B2;
    % Computing B3 with operation: mult
    B3=B2 * Bb3;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.4999999999999432 + 0.0im;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=1.0 + 0.0im;
    T2k5= coeff1*T2k3+coeff2*B3;
    output=T2k5;
end

