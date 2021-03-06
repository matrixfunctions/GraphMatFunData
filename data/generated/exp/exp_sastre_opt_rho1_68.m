function output = exp_sastre_opt_rho1_68(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1 = -0.3371652952117645;
    coeff2 = 0.4456299027272243;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = -0.2940691417403946;
    coeff2 = 0.4830197021196582;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = 0.0007801504959769123;
    coeff2 = 0.036934048204913826;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 0.008399213058466935;
    coeff2 = 0.0017935941807155679;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -0.2940691417403946;
    coeff2 = 0.4830197021196582;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.08435039764895878;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.33589987345926525;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0005218334570419188;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = 0.968291126022028;
    coeff2 = 0.5165851751813124;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.47391353450112705;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = 0.11409149134408365;
    coeff2 = 0.68729265604707;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.37829567639304673;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 0.1382844980970843;
    coeff2 = 0.13000362502263652;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.06285653463395623;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1 = -0.055527011448072344;
    coeff2 = -0.019592818624400106;
    T2k2 = coeff1*I + coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.01516878804519063;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 0.968291126022028;
    coeff2 = 0.5165851751813124;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.47391353450112705;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 0.04182430699533517;
    coeff2 = -0.19033210762698088;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9365591932640623;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9668380398527914;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 7.480755260549251e-5;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.11560763975481574;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.033903885370555546;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.36476063896683886;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.160367463472582;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.002839711497134736;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.014357966768578765;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.018703010057794192;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.11560763975481574;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6395504919932988;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.160367463472582;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6818767988800816;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6818767988800816;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.009741431366481671;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0239727668224385;
    T2k8 = coeff1*T2k6 + coeff2*B6;
    output = T2k8;
end

