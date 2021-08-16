function output = exp_sastre_m5_opt_rho1_9(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1 = -0.35706341164511884;
    coeff2 = 0.41479581112660663;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = -0.29597130380221776;
    coeff2 = 0.4830650480818557;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -0.005411040008503428;
    coeff2 = 0.0479465268958071;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 0.014197999026813905;
    coeff2 = 0.0023989722622025037;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -0.29597130380221776;
    coeff2 = 0.4830650480818557;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.08232725511150228;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.3539095385982056;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0004504909091063826;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = 1.0985965640823765;
    coeff2 = 0.45359605154012905;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.3807951389639113;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = 0.10677883331398988;
    coeff2 = 0.6584173995757578;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.4164266931927851;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 0.12273251114319797;
    coeff2 = 0.1799587581454721;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.043513380168846295;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1 = -0.050949343525589726;
    coeff2 = -0.019581395624785863;
    T2k2 = coeff1*I + coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.016433672152666983;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 0.9068877097254958;
    coeff2 = 0.5793395460985866;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5479499415010913;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 0.052794118729278604;
    coeff2 = -0.22109848126655623;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9222178419369582;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9600894604220988;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00011985947592822976;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.13249960107719003;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03716863529834358;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.2855172295604338;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.153379110061045;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0018851674132650941;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0009698855603010127;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.02122483727644592;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.09228091229564235;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5818565692388099;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.179587866680269;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6801751580137498;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6370429665254512;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.014316598391468338;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9728079294146128;
    T2k8 = coeff1*T2k6 + coeff2*B6;
    output = T2k8;
end

