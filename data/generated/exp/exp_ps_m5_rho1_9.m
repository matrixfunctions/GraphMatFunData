function output = exp_ps_m5_rho1_9(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Bb5_2 T2k2 Bb6_2 B2 T2k3 Bb6_3 Bb5_3 B3 Bb6_4 T2k4 Bb5_4 B4 Bb5 B5 Bb6 B6 T2k8
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 2.48015873015873e-5;
    coeff2 = 2.7557319223985893e-6;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 0.041666666666666664;
    coeff2 = 0.008333333333333333;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.001388888888888889;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.755731922398589e-7;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = B2 * A;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0001984126984126984;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.16666666666666666;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.505210838544172e-8;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = B3 * A;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.08767569878681e-9;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = B4 * Bb5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Bb6 = coeff1*Bb6_4 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = B4 * Bb6;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k8 = coeff1*T2k4 + coeff2*B6;
    output = T2k8;
end

