function output = exp_bbcs_m3(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Bb4_2 T2k2 B2 Bb3 T2k3 Bb4_3 B3 Ba4 Bb4 B4 T2k6
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 0.0 + 1i*0.5496085391143601;
    coeff2 = 0.1620095284677366 + 0.0im;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 1.0 + 0.0im;
    coeff2 = -0.0 + 1i*-0.9999999999999923;
    T2k2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 0.10775 + 0.0im;
    coeff2 = -0.0 + 1i*-0.026939068735988708;
    Bb3 = coeff1*A + coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0 + 0.0im;
    coeff2 = -0.13549409636220702 + 0.0im;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0 + 0.0im;
    coeff2 = -0.0 + 1i*-0.014179818052118045;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = B2 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 0.0 + 1i*0.6632100444166243;
    coeff2 = 1.0 + 0.0im;
    Ba4 = coeff1*B2 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0 + 0.0im;
    coeff2 = -0.034159539168921116 + 0.0im;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0 + 0.0im;
    coeff2 = 1.0 + 0.0im;
    T2k6 = coeff1*T2k3 + coeff2*B4;
    output = T2k6;
end

