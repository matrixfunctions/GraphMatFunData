function output = exp_mono_m6_opt_rho2_7(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1 = -0.039878213296149906;
    coeff2 = 0.3191312115175403;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = -0.000118740658609087;
    coeff2 = -0.00020274697859323616;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 0.9837654928593774;
    coeff2 = 0.8167945494192191;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = 0.1980540938521463;
    coeff2 = 0.9168156264182638;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -0.012054524614470144;
    coeff2 = -0.015473430203696856;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = -0.011642655100621168;
    coeff2 = 1.0007182791901643;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = 0.1980540938521463;
    coeff2 = 0.9168156264182638;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.3230094878143697;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9305465620354068;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.004184759286380629;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.41828764759110026;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.05844498049251375;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -0.001434162774737488;
    coeff2 = -0.005882107697567988;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.004345083585047115;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 0.8274052334396728;
    coeff2 = 0.6467470861357121;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.052798886182833585;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = -4.49699371062923e-6;
    coeff2 = 4.300156065032907e-6;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.354896741415178e-5;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = -0.10696304090753009;
    coeff2 = 0.9888245872897673;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.10240241498554757;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 0.30406437703941375;
    coeff2 = 0.9357055778901153;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.031045948902277486;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9639291036056933;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00182483121117615;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00845811634632641;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.013669299358486304;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.8211621256077539;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0002567501403605987;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1592701070085137;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03329710767309;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0030548445338270413;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.39542321609363e-5;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00850052220135662;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.308226007116608e-7;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.13046443012976514;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.8480144820381275;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.8140979919427691;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0010349274617070035;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = -8.044483063015104e-8;
    coeff2 = 0.9965080883913857;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.11832377144773065;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.011665314654668413;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.00613552979940033;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0005336901003739011;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.06890626823439797;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.726903525987751e-5;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.12050842503745675;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9947846525992833;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.4795445648642135e-8;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.1578020415467386e-5;
    T2k9 = coeff1*T2k7 + coeff2*B7;
    output = T2k9;
end

