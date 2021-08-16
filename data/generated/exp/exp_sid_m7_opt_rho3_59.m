function output = exp_sid_m7_opt_rho3_59(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 0.0005326111625306977;
    coeff2 = 1.0012299933973516;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = -1.5738212245123645e-8;
    coeff2 = 3.740920280371873e-8;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 0.9999420628667928;
    coeff2 = 0.9990039337658048;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = 0.00022214082957652806;
    coeff2 = 5.947941137388186;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = -0.00016401351885017156;
    coeff2 = 0.0005433816370422478;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = -0.0017737608141714808;
    coeff2 = 1.005235674329508;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = 9.544368792763562e-5;
    coeff2 = 0.9999194856569542;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = -1.1883273082175028e-5;
    coeff2 = 3.080294633483493e-5;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1 = 0.00028017073414700555;
    coeff2 = -0.0013783167325222296;
    Bb8_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -0.0017737608141714808;
    coeff2 = 1.005235674329508;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0012270195002069;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.007816923517808423;
    Bb8_3 = coeff1*Bb8_2 + coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0036194298213211165;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.006925033919227824;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.4220821708139979;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.024845741420535115;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.00033449839479958844;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.267964654623867e-7;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0011284750735624438;
    Bb8_4 = coeff1*Bb8_3 + coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5079542977446111;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.08249596102472202;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 9.623360783744299e-5;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9998936828701804;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 8.405278388325319e-5;
    Bb8_5 = coeff1*Bb8_4 + coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.005411384044110833;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -3.46635995191657e-6;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.708627877122971e-6;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0011101347889166462;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -6.0508143407070435e-6;
    coeff2 = 1.0000102174714505;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0028500766374399016;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.00013226806180853073;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0208196110243714e-5;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = -4.786848660402673e-7;
    coeff2 = 8.582321314936397e-8;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.9167378426502e-6;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -5.991737523333378e-6;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.999999999087183;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 7.400879751137793e-5;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999999999895876;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.168781522592016e-5;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.3621609668499039e-5;
    Bb8_6 = coeff1*Bb8_5 + coeff2*B5;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = -6.288899284895006e-6;
    coeff2 = 1.1355716740441329e-6;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.567958735815622e-8;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.487015589947048e-9;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0836272945889375e-10;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.1368496406363211e-11;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.003633457505924;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -4.7284436864949066e-5;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9977252838919913;
    Bb8_7 = coeff1*Bb8_6 + coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1 = 0.2015451143023047;
    coeff2 = 11.106844258571114;
    Ba8_2 = coeff1*I + coeff2*A;
    % Computing Ba8_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.9912660184614275;
    Ba8_3 = coeff1*Ba8_2 + coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.2891134043244349;
    Ba8_4 = coeff1*Ba8_3 + coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.02649378029184827;
    Ba8_5 = coeff1*Ba8_4 + coeff2*B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0027421911866186664;
    Ba8_6 = coeff1*Ba8_5 + coeff2*B5;
    % Computing Ba8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 61.7597688836338;
    Ba8_7 = coeff1*Ba8_6 + coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = -0.00014732216151790324;
    coeff2 = -0.004469043351396654;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0326061354999038;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0028514008243872304;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0001864934813936358;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.852437481436439e-5;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9846670594297366;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9846397447218478;
    Ba8 = coeff1*Ba8_7 + coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.4108017105321994e-5;
    Bb8 = coeff1*Bb8_7 + coeff2*B7;
    % Computing B8 with operation: mult
    B8 = Ba8 * Bb8;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0023105536908983168;
    T2k8 = coeff1*T2k7 + coeff2*B7;
    % Computing T2k10 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9977227168431884;
    T2k10 = coeff1*T2k8 + coeff2*B8;
    output = T2k10;
end

