function output = exp_ps_m9(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: B_3_1 B_0_1 B_4_1 B_1_1 B_2_1 A2 B_4_2 B_0_2 B_2_2 B_3_2 B_1_2 A3 B_2_3 B_3_3 B_0_3 B_4_3 B_1_3 A4 B_3_4 B_4_4 B_2_4 B_1_4 B_0_4 A5 B_4_5 B_1_5 B_3_5 B_0_5 B_2_5 A6 P4 C3 P3 C2 P2 C1 P1 C0 P0
    % Computing B_3_1 with operation: lincomb
    coeff1 = 1.5619206968586225e-16;
    coeff2 = 8.22063524662433e-18;
    B_3_1 = coeff1*I + coeff2*A;
    % Computing B_0_1 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    B_0_1 = coeff1*I + coeff2*A;
    % Computing B_4_1 with operation: lincomb
    coeff1 = 1.6117375710961184e-24;
    coeff2 = 6.446950284384474e-26;
    B_4_1 = coeff1*I + coeff2*A;
    % Computing B_1_1 with operation: lincomb
    coeff1 = 0.001388888888888889;
    coeff2 = 0.0001984126984126984;
    B_1_1 = coeff1*I + coeff2*A;
    % Computing B_2_1 with operation: lincomb
    coeff1 = 2.08767569878681e-9;
    coeff2 = 1.6059043836821613e-10;
    B_2_1 = coeff1*I + coeff2*A;
    % Computing A2 with operation: mult
    A2 = A * A;
    % Computing B_4_2 with operation: lincomb
    coeff1 = 2.4795962632247976e-27;
    coeff2 = 1.0;
    B_4_2 = coeff1*A2 + coeff2*B_4_1;
    % Computing B_0_2 with operation: lincomb
    coeff1 = 0.5;
    coeff2 = 1.0;
    B_0_2 = coeff1*A2 + coeff2*B_0_1;
    % Computing B_2_2 with operation: lincomb
    coeff1 = 1.1470745597729725e-11;
    coeff2 = 1.0;
    B_2_2 = coeff1*A2 + coeff2*B_2_1;
    % Computing B_3_2 with operation: lincomb
    coeff1 = 4.110317623312165e-19;
    coeff2 = 1.0;
    B_3_2 = coeff1*A2 + coeff2*B_3_1;
    % Computing B_1_2 with operation: lincomb
    coeff1 = 2.48015873015873e-5;
    coeff2 = 1.0;
    B_1_2 = coeff1*A2 + coeff2*B_1_1;
    % Computing A3 with operation: mult
    A3 = A * A2;
    % Computing B_2_3 with operation: lincomb
    coeff1 = 7.647163731819816e-13;
    coeff2 = 1.0;
    B_2_3 = coeff1*A3 + coeff2*B_2_2;
    % Computing B_3_3 with operation: lincomb
    coeff1 = 1.9572941063391263e-20;
    coeff2 = 1.0;
    B_3_3 = coeff1*A3 + coeff2*B_3_2;
    % Computing B_0_3 with operation: lincomb
    coeff1 = 0.16666666666666666;
    coeff2 = 1.0;
    B_0_3 = coeff1*A3 + coeff2*B_0_2;
    % Computing B_4_3 with operation: lincomb
    coeff1 = 9.183689863795546e-29;
    coeff2 = 1.0;
    B_4_3 = coeff1*A3 + coeff2*B_4_2;
    % Computing B_1_3 with operation: lincomb
    coeff1 = 2.7557319223985893e-6;
    coeff2 = 1.0;
    B_1_3 = coeff1*A3 + coeff2*B_1_2;
    % Computing A4 with operation: mult
    A4 = A * A3;
    % Computing B_3_4 with operation: lincomb
    coeff1 = 8.896791392450574e-22;
    coeff2 = 1.0;
    B_3_4 = coeff1*A4 + coeff2*B_3_3;
    % Computing B_4_4 with operation: lincomb
    coeff1 = 3.279889237069838e-30;
    coeff2 = 1.0;
    B_4_4 = coeff1*A4 + coeff2*B_4_3;
    % Computing B_2_4 with operation: lincomb
    coeff1 = 4.779477332387385e-14;
    coeff2 = 1.0;
    B_2_4 = coeff1*A4 + coeff2*B_2_3;
    % Computing B_1_4 with operation: lincomb
    coeff1 = 2.755731922398589e-7;
    coeff2 = 1.0;
    B_1_4 = coeff1*A4 + coeff2*B_1_3;
    % Computing B_0_4 with operation: lincomb
    coeff1 = 0.041666666666666664;
    coeff2 = 1.0;
    B_0_4 = coeff1*A4 + coeff2*B_0_3;
    % Computing A5 with operation: mult
    A5 = A * A4;
    % Computing B_4_5 with operation: lincomb
    coeff1 = 1.1309962886447716e-31;
    coeff2 = 1.0;
    B_4_5 = coeff1*A5 + coeff2*B_4_4;
    % Computing B_1_5 with operation: lincomb
    coeff1 = 2.505210838544172e-8;
    coeff2 = 1.0;
    B_1_5 = coeff1*A5 + coeff2*B_1_4;
    % Computing B_3_5 with operation: lincomb
    coeff1 = 3.868170170630684e-23;
    coeff2 = 1.0;
    B_3_5 = coeff1*A5 + coeff2*B_3_4;
    % Computing B_0_5 with operation: lincomb
    coeff1 = 0.008333333333333333;
    coeff2 = 1.0;
    B_0_5 = coeff1*A5 + coeff2*B_0_4;
    % Computing B_2_5 with operation: lincomb
    coeff1 = 2.8114572543455206e-15;
    coeff2 = 1.0;
    B_2_5 = coeff1*A5 + coeff2*B_2_4;
    % Computing A6 with operation: mult
    A6 = A * A5;
    % Computing P4 with operation: lincomb
    coeff1 = 3.7699876288159054e-33;
    coeff2 = 1.0;
    P4 = coeff1*A6 + coeff2*B_4_5;
    % Computing C3 with operation: mult
    C3 = P4 * A6;
    % Computing P3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    P3 = coeff1*C3 + coeff2*B_3_5;
    % Computing C2 with operation: mult
    C2 = P3 * A6;
    % Computing P2 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    P2 = coeff1*C2 + coeff2*B_2_5;
    % Computing C1 with operation: mult
    C1 = P2 * A6;
    % Computing P1 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    P1 = coeff1*C1 + coeff2*B_1_5;
    % Computing C0 with operation: mult
    C0 = P1 * A6;
    % Computing P0 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    P0 = coeff1*C0 + coeff2*B_0_5;
    output = P0;
end

