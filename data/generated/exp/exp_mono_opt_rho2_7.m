function output = exp_mono_opt_rho2_7(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1 = -3.857235106310176;
    coeff2 = 2.559191063273799;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = -0.027531215857697383;
    coeff2 = 1.4691934416553218;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 9.22780659770276;
    coeff2 = -7.818752581614379;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = -0.8990318926117189;
    coeff2 = 1.0529942676288555;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -34.51931508412978;
    coeff2 = -32.102642943446725;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 22.563420733949638;
    coeff2 = 41.22336752504774;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -0.8990318926117189;
    coeff2 = 1.0529942676288555;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.7030375734689573;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.2675261675885479;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.327421400455049;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -4.7517129085184555;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 8.606161607327088;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -1.6333763880809702;
    coeff2 = -12.425428631305296;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -7.020656636475915;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 2.576203242554731;
    coeff2 = 10.289018418102025;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.34028419625059303;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = 0.20156975628485335;
    coeff2 = -0.11857297824908863;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.13933557327081972;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = -1.0983539425730369;
    coeff2 = 4.121404800758164;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.10671574515634609;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 2.6127319091880588;
    coeff2 = -3.5281474819181633;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -4.966605691968948;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.5107167490732225;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 18.485870531498684;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5170030863801119;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 15.340351060111168;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.961642608409374;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -2.484455273478521;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 29.297986883148763;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0069810802300877255;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 26.83688350147083;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -12.951270741270744;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0011903433805721666;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -5.578338573600123;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.5545952756429124;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.47775953255938497;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -4.727574599810485;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.1317145971276296e-9;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = 1.0897581376826212e-7;
    coeff2 = -1.555574375395291;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.8159169191071185;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.13432676772501073;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.008988745081115584;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.5651632231929153;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.6798469796207142;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.3062199878200116;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.5814827606241575;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.5449585814160118;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.022514410581526905;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.3350038182653414e-45;
    T2k9 = coeff1*T2k7 + coeff2*B7;
    output = T2k9;
end

