function output = exp_sastre_m8_opt_rho13_5(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Ba4_2 Ba9_2 Bb9_2 Bb2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba7_3 Ba9_3 Ba3 B3 Ba7_4 Ba4 Bb4_3 Bb9_3 Bb8_4 Bb4 B4 Bb8_5 Ba6_3 Ba9_4 Ba9_5 T2k3 T2k4 Bb9_4 Bb9_5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Bb9_6 Ba6 Ba9_6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Bb9_7 Ba7 Ba9_7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba9_8 Bb9_8 Ba8 Bb8 T2k8 B8 Ba9 T2k9 Bb9 B9 T2k11
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 0.004105729412589247;
    coeff2 = -0.006768382232040512;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = 0.011539085894536814;
    coeff2 = 0.0935630221089867;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 0.08072764899462631;
    coeff2 = 0.02311006099515489;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = 3.0826131069323353;
    coeff2 = 0.2781411205676187;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = -0.0028647512623116655;
    coeff2 = 0.0009954035208072304;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = 0.027019231244596563;
    coeff2 = 0.06170763799314927;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Ba9_2 with operation: lincomb
    coeff1 = 0.3035081839334926;
    coeff2 = 0.042128860424959606;
    Ba9_2 = coeff1*I + coeff2*A;
    % Computing Bb9_2 with operation: lincomb
    coeff1 = -0.2738445106644693;
    coeff2 = -0.033002263370480726;
    Bb9_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = -0.11490308438229017;
    coeff2 = 0.15730697013710307;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = -0.013620620860785195;
    coeff2 = 0.07638110624385147;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1 = 0.9222548073365477;
    coeff2 = 0.10238368416498037;
    Bb8_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -0.11490308438229017;
    coeff2 = 0.15730697013710307;
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
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.15593200742501928;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Ba9_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.005532874783353146;
    Ba9_3 = coeff1*Ba9_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999471251526917;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 7.925112972476668;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0067849094795742;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.14911649238384225;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Bb9_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.004912635491924672;
    Bb9_3 = coeff1*Bb9_2 + coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.030058902112840666;
    Bb8_4 = coeff1*Bb8_3 + coeff2*B3;
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
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.06407916439324149;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Ba9_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.05011906670785025;
    Ba9_4 = coeff1*Ba9_3 + coeff2*B3;
    % Computing Ba9_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03702751170872807;
    Ba9_5 = coeff1*Ba9_4 + coeff2*B4;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.05527883142327172;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.5378493201072977e-5;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb9_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.04449492500250622;
    Bb9_4 = coeff1*Bb9_3 + coeff2*B3;
    % Computing Bb9_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03287724382394223;
    Bb9_5 = coeff1*Bb9_4 + coeff2*B4;
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
    coeff2 = -0.13204423967508136;
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
    coeff2 = 0.001517105603937254;
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
    coeff2 = 1.292606058383822e-5;
    T2k6 = coeff1*T2k4 + coeff2*B5;
    % Computing Bb9_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.05879027030590609;
    Bb9_6 = coeff1*Bb9_5 + coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0008083736633122;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Ba9_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.06620554678336489;
    Ba9_6 = coeff1*Ba9_5 + coeff2*B5;
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
    coeff2 = 0.016824149148287185;
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
    % Computing Bb9_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.13983408425497476;
    Bb9_7 = coeff1*Bb9_6 + coeff2*B6;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.16056529744528691;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing Ba9_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.15748603557212862;
    Ba9_7 = coeff1*Ba9_6 + coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 6.578027159862427e-7;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.05247830183077817;
    Bb8_7 = coeff1*Bb8_6 + coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1 = 0.9222548073365477;
    coeff2 = 0.10238368416498037;
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
    coeff2 = 0.05228210309031749;
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
    % Computing Ba9_8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.03947847759964156;
    Ba9_8 = coeff1*Ba9_7 + coeff2*B7;
    % Computing Bb9_8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.03505336576885022;
    Bb9_8 = coeff1*Bb9_7 + coeff2*B7;
    % Computing Ba8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9972225209602774;
    Ba8 = coeff1*Ba8_7 + coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9972225209602774;
    Bb8 = coeff1*Bb8_7 + coeff2*B7;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.3533497256574092e-7;
    T2k8 = coeff1*T2k7 + coeff2*B7;
    % Computing B8 with operation: mult
    B8 = Ba8 * Bb8;
    % Computing Ba9 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.058382881818128;
    Ba9 = coeff1*Ba9_8 + coeff2*B8;
    % Computing T2k9 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.9207658227240783e-6;
    T2k9 = coeff1*T2k8 + coeff2*B8;
    % Computing Bb9 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9397507089701836;
    Bb9 = coeff1*Bb9_8 + coeff2*B8;
    % Computing B9 with operation: mult
    B9 = Ba9 * Bb9;
    % Computing T2k11 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9800057467529937;
    T2k11 = coeff1*T2k9 + coeff2*B9;
    output = T2k11;
end

