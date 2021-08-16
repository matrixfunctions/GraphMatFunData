function output = exp_mono_m7_opt_rho3_59(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10
    % Computing Ba3_2 with operation: lincomb
    coeff1 = -8.973017805624113e-5;
    coeff2 = 0.00012005201731132445;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = -3.328521860384963e-7;
    coeff2 = -2.639951734327968e-6;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 1.00000000175563;
    coeff2 = 1.0000017846980391;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = -2.0909266321867902e-8;
    coeff2 = -2.7287925158899476e-7;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 0.00012338976185732696;
    coeff2 = 0.9996312949109102;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = 1.3134712061518662e-5;
    coeff2 = 1.002274766505785;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -1.8416000490973626e-5;
    coeff2 = -8.104101755291844e-5;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 4.02314359582651e-6;
    coeff2 = 0.9993459669804694;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1 = -1.0645153851732277e-9;
    coeff2 = 0.9987013218201104;
    Bb8_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = 1.3134712061518662e-5;
    coeff2 = 1.002274766505785;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.018858289778215145;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.052949897775355415;
    Bb8_3 = coeff1*Bb8_2 + coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.489518987562822e-5;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9997367268186822;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -8.505783380079528e-8;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -2.4699349369528174e-7;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9998500554706558;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.1160650661114392e-6;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.003442617642989778;
    Bb8_4 = coeff1*Bb8_3 + coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.49769531969924596;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.16600417761584907;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03609843676884375;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0010584605114307261;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00011084212403663171;
    Bb8_5 = coeff1*Bb8_4 + coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03838601189631294;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.185994260757166e-6;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 7.151219728856269e-5;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 7.547096805172991e-6;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -2.5675820144906785e-6;
    coeff2 = -1.618163016765615e-5;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -3.4221784684677797e-6;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.4268115376905545e-5;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999837106332966;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 4.578521649969897e-5;
    coeff2 = 1.0001351072704798;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.013461797467855375;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.002324571658090968;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0002315295921500602;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.006185920919246356;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999995656134203;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.00011038627353830762;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -3.8270257669926835e-7;
    Bb8_6 = coeff1*Bb8_5 + coeff2*B5;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = -0.00011293613093022973;
    coeff2 = 0.9990923104443914;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.047438292970867935;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0033538018355920127;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00011542577452508812;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -2.150293392583585e-7;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999999715919681;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0010600598008071952;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.4340108650239077e-8;
    Bb8_7 = coeff1*Bb8_6 + coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1 = 0.40308943704974315;
    coeff2 = -1.1533961194636121e-8;
    Ba8_2 = coeff1*I + coeff2*A;
    % Computing Ba8_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -5.03389309270026e-9;
    Ba8_3 = coeff1*Ba8_2 + coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.635242146052482e-8;
    Ba8_4 = coeff1*Ba8_3 + coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.345288563528515e-7;
    Ba8_5 = coeff1*Ba8_4 + coeff2*B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 8.355043856393377e-7;
    Ba8_6 = coeff1*Ba8_5 + coeff2*B5;
    % Computing Ba8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.00012492568018016112;
    Ba8_7 = coeff1*Ba8_6 + coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = -0.00012498199105449398;
    coeff2 = 0.9989757606293848;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.04754175140634288;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.002515599835204672;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 8.584302495182615e-5;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.619813578008144e-7;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 7.644984585045822e-9;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9999999932294968;
    Ba8 = coeff1*Ba8_7 + coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 9.818180030024404e-11;
    Bb8 = coeff1*Bb8_7 + coeff2*B7;
    % Computing B8 with operation: mult
    B8 = Ba8 * Bb8;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 8.103038123336939e-5;
    T2k8 = coeff1*T2k7 + coeff2*B7;
    % Computing T2k10 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.7925485378650857e-6;
    T2k10 = coeff1*T2k8 + coeff2*B8;
    output = T2k10;
end

