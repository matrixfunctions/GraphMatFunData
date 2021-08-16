function output = exp_sastre_m7_opt_rho6_4(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 0.004105729412589247;
    coeff2 = -0.013536764464081023;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = 0.011539085894536814;
    coeff2 = 0.1871260442179734;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = -0.002287189805563291;
    coeff2 = 0.004627231318175632;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = 3.0826131069323353;
    coeff2 = 0.5562822411352374;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = -0.0028647512623116655;
    coeff2 = 0.001990807041614461;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = -0.11490308438229017;
    coeff2 = 0.31461394027420614;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = 0.027019231244596563;
    coeff2 = 0.12341527598629853;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = -0.013620620860785195;
    coeff2 = 0.15276221248770294;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1 = 0.9222548073365477;
    coeff2 = 0.20476736832996073;
    Bb8_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -0.11490308438229017;
    coeff2 = 0.31461394027420614;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.001211033078915847;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1032554965257702;
    Bb8_3 = coeff1*Bb8_2 + coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0018954361862412274;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999471251526917;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.15593200742501928;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 7.925112972476668;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0067849094795742;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.06407916439324149;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.030058902112840666;
    Bb8_4 = coeff1*Bb8_3 + coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.005161155569679917;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.046748826502193745;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.14911649238384225;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9051158545637227;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.11384698074301966;
    Bb8_5 = coeff1*Bb8_4 + coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.034540156044914797;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.285986613635287;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.044725605639523214;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.001825018602058;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -0.0032451095696861268;
    coeff2 = -0.2640884793501627;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9186197117187667;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.42211811124544796;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.01259382817947474;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 0.011274427044379665;
    coeff2 = 0.003034211207874508;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.001158435856997462;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0047467287622741695;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -7.137183226929846e-5;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.061761015137031604;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0008083736633122;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.10421074361186333;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.011734571936531983;
    Bb8_6 = coeff1*Bb8_5 + coeff2*B5;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 0.02459649940876415;
    coeff2 = 0.03364829829657437;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03458106185410125;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.047464943469750574;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.11471274073059659;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0130277788963542;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.16056529744528691;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.14690678556135978;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.05247830183077817;
    Bb8_7 = coeff1*Bb8_6 + coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1 = 0.9222548073365477;
    coeff2 = 0.20476736832996073;
    Ba8_2 = coeff1*I + coeff2*A;
    % Computing Ba8_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1032554965257702;
    Ba8_3 = coeff1*Ba8_2 + coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.030058902112840666;
    Ba8_4 = coeff1*Ba8_3 + coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.11384698074301966;
    Ba8_5 = coeff1*Ba8_4 + coeff2*B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.011734571936531983;
    Ba8_6 = coeff1*Ba8_5 + coeff2*B5;
    % Computing Ba8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.05247830183077817;
    Ba8_7 = coeff1*Ba8_6 + coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = 0.025970161938770917;
    coeff2 = 0.10456420618063499;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.10608150003582509;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.049985371773829376;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.28766922644726883;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.936009578758007;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0130307509989027;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9972225209602774;
    Ba8 = coeff1*Ba8_7 + coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9972225209602774;
    Bb8 = coeff1*Bb8_7 + coeff2*B7;
    % Computing B8 with operation: mult
    B8 = Ba8 * Bb8;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.03682640896812769;
    T2k8 = coeff1*T2k7 + coeff2*B7;
    % Computing T2k10 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9872838791818296;
    T2k10 = coeff1*T2k8 + coeff2*B8;
    output = T2k10;
end

