function output = exp_ps_m6_opt_rho2_7(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 0.19178054915267911;
    coeff2 = 0.7596053740559234;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = 0.0006283647624375745;
    coeff2 = -0.0856504866246125;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 0.6845260319725605;
    coeff2 = 0.1283513925026888;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = 0.729254790966769;
    coeff2 = 0.36467163874605407;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -0.0050940656098555845;
    coeff2 = 0.3344727118113129;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 1.2113351965247428;
    coeff2 = 0.08511915910470753;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = 0.7280078546464023;
    coeff2 = 0.3646378676083862;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9346441756079398;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.722565841223338;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.06273050414017863;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.021508466841232733;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0039775493520180115;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = 1.3379066745135616e-7;
    coeff2 = -0.002462022732959875;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03688690867443445;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 2.3798863906382593e-5;
    coeff2 = 1.0837308004804432e-6;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.2489586532561667e-6;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = 0.9919381605205452;
    coeff2 = 0.26220671903598886;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.26972239418076976;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 0.009818461198297107;
    coeff2 = 0.0010390323767811615;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.000350687159799907;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 1.0908348608940441;
    coeff2 = 0.26673427641075104;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1599439900856757;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.789504699569824;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.006605177298124698;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 8.23863768936227e-5;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.20340771723284812;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.986721033920639;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.618679223329589e-5;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.055068936428514716;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.2914936642490679;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.0091806253298255e-9;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.380604588051657e-8;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0432353938528486e-5;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999998069262623;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.8231044911539176;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.954079109272231;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.689858409164275e-5;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00037551432005694267;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = -0.05660472123179401;
    coeff2 = 0.021982855815076786;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.032263911106579925;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.021016177521655365;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5440682102221873;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.021400262174087047;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.01874072744531;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0006732366416760824;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.3499625112901884e-5;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.05804569963149164;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0044248526362387;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0045369223085259;
    T2k9 = coeff1*T2k7 + coeff2*B7;
    output = T2k9;
end

