function output = exp_ps_m7_rho6_0(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Bb8_2 T2k2 Bb6_2 Bb7_2 B2 Bb7_3 T2k3 Bb6_3 Bb8_3 B3 Bb6_4 Bb7_4 T2k4 Bb8_4 B4 Bb8_5 T2k5 Bb6_5 Bb7_5 B5 Bb6 B6 Bb7 B7 Bb8 Ba8_6 B8 T2k10
    % Computing Bb8_2 with operation: lincomb
    coeff1 = 0.008333333333333333;
    coeff2 = 0.001388888888888889;
    Bb8_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 7.647163731819816e-13;
    coeff2 = 4.779477332387385e-14;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = 2.755731922398589e-7;
    coeff2 = 2.505210838544172e-8;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.08767569878681e-9;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.8114572543455206e-15;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0001984126984126984;
    Bb8_3 = coeff1*Bb8_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = B2 * A;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.5619206968586225e-16;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.6059043836821613e-10;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.16666666666666666;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.48015873015873e-5;
    Bb8_4 = coeff1*Bb8_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = B3 * A;
    % Computing Bb8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.7557319223985893e-6;
    Bb8_5 = coeff1*Bb8_4 + coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.041666666666666664;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 8.22063524662433e-18;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.1470745597729725e-11;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = B4 * A;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.110317623312165e-19;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = B5 * Bb6;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Bb7 = coeff1*Bb7_5 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = B5 * Bb7;
    % Computing Bb8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Bb8 = coeff1*Bb8_5 + coeff2*B7;
    % Computing Ba8_6 with operation: lincomb
    coeff1 = 0.40308943770549005;
    coeff2 = 1.0;
    Ba8_6 = coeff1*I + coeff2*B5;
    % Computing B8 with operation: mult
    B8 = Ba8_6 * Bb8;
    % Computing T2k10 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k10 = coeff1*T2k5 + coeff2*B8;
    output = T2k10;
end

