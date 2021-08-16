function output = exp_sid_m5_opt_rho1_68(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 4.5138550333531274e-8;
    coeff2 = 0.9999996252970185;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = -5.497501731236461e-7;
    coeff2 = 0.9999997231244894;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -9.77343257577906e-10;
    coeff2 = -1.0901600119737703e-9;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = -1.5828863347334606e-7;
    coeff2 = 5.3740768955983805e-5;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -5.497501731236461e-7;
    coeff2 = 0.9999997231244894;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 9.489090720230229e-7;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -7.203565312107398e-12;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.501624603111096e-6;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = 3.3031893066419052e-6;
    coeff2 = 0.31120539066317016;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.08528254475386922;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = 7.960623161921035e-6;
    coeff2 = 0.9418651592948695;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.06974385503852226;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = -1.5696506919081925e-6;
    coeff2 = 3.82192171574818e-6;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.007543497631342593;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1 = 1.0000000000279976;
    coeff2 = 1.0000016922101196;
    T2k2 = coeff1*I + coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.2863259306351655;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 3.7295572226790776e-6;
    coeff2 = 0.6865711102982829;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.13922314499608707;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = -2.9371772095919626e-7;
    coeff2 = 4.5138353687845936e-8;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999996252973244;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.999999999999673;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.161609126932939e-6;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.031514091457595775;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0020055991068182832;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.999989585781789;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 11.17362496085854;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.08775648076729538;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1899551457347289;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0028530586220509276;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.029243268370171514;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999932366972974;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.8297817176379367;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999873285367696;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999948843658316;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.2333735694151144;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999971446794519;
    T2k8 = coeff1*T2k6 + coeff2*B6;
    output = T2k8;
end

