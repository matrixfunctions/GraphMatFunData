function output = exp_ps_m7_opt_rho6_4(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 0.09654809779701808;
    coeff2 = 0.9314722195252554;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = -2.47688825760804e-9;
    coeff2 = 0.0004847383569867066;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 0.7602626115113443;
    coeff2 = 0.2868135868803526;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = 5.241383560289876e-5;
    coeff2 = 0.014487073789250026;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 0.9062397883450928;
    coeff2 = 0.11339252602614255;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = 0.42031771622120523;
    coeff2 = 0.309252105381124;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -0.06178357797354583;
    coeff2 = 0.13477081084608053;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 0.7648819021317048;
    coeff2 = 0.3110254108862011;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1 = 0.09171588796170577;
    coeff2 = 0.01956469387438758;
    Bb8_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = 0.42031771622120523;
    coeff2 = 0.309252105381124;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.12302744324936216;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.01590736865892108;
    Bb8_3 = coeff1*Bb8_2 + coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.7329752440067754;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.19731022777550086;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.01054251273835961;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.12590635459023522;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5453304515439815;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.008759615364757254;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.010821193755227018;
    Bb8_4 = coeff1*Bb8_3 + coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.36416792067469345;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1470803112875071;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.13872927228389606;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.031072070853003377;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.002975380215829149;
    Bb8_5 = coeff1*Bb8_4 + coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.301249740245868;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.024971360914667298;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.20141416233600964;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.35076815913299725;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -0.0037226968956481066;
    coeff2 = 0.3027942237668864;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.16951280499357504;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.7173684223758208;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6151562256858113;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 0.4984122280039027;
    coeff2 = 0.6168934079023138;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.32464232547450916;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.11022586483948126;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00562561032811143;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.05977870243645669;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9838034227912541;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9784008568588043;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.005577510585906947;
    Bb8_6 = coeff1*Bb8_5 + coeff2*B5;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = -7.332920860083383e-8;
    coeff2 = 7.255106883862399e-10;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.599720660809858e-11;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.386947932212966e-11;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.696432230206371e-12;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.7604767541951467e-14;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.756683002701541e-7;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.4088571170380866e-7;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -7.0435154500906185e-6;
    Bb8_7 = coeff1*Bb8_6 + coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1 = 0.9900217594371604;
    coeff2 = -0.0426544708742059;
    Ba8_2 = coeff1*I + coeff2*A;
    % Computing Ba8_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0895038182070465;
    Ba8_3 = coeff1*Ba8_2 + coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.4332397850086141;
    Ba8_4 = coeff1*Ba8_3 + coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6230790774477177;
    Ba8_5 = coeff1*Ba8_4 + coeff2*B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.836798518267107;
    Ba8_6 = coeff1*Ba8_5 + coeff2*B5;
    % Computing Ba8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 6.715240076479665e-7;
    Ba8_7 = coeff1*Ba8_6 + coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = -0.004732216406789163;
    coeff2 = 2.4676441956771072e-5;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.805471807800756e-5;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 6.965400723989983e-6;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.4212591399083807e-6;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.2164967481538e-7;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0000000002008125;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03410244422083362;
    Ba8 = coeff1*Ba8_7 + coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0199299719677057;
    Bb8 = coeff1*Bb8_7 + coeff2*B7;
    % Computing B8 with operation: mult
    B8 = Ba8 * Bb8;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00148688702927498;
    T2k8 = coeff1*T2k7 + coeff2*B7;
    % Computing T2k10 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.021358342405418;
    T2k10 = coeff1*T2k8 + coeff2*B8;
    output = T2k10;
end

