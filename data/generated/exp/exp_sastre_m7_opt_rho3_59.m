function output=exp_sastre_m7_opt_rho3_59(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10
    % Computing Ba3_2 with operation: lincomb
    coeff1=0.00203408685542654;
    coeff2=-0.0024277753235612716;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=0.01563610996898412;
    coeff2=0.18417388850684976;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=5.747397676949589e-7;
    coeff2=2.761844663337101e-7;
    T2k2= coeff1*I+coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1=3.1043462419448544;
    coeff2=0.3764280070972597;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1=-0.006811245717889387;
    coeff2=0.0009029116501382352;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.0014416692554755693;
    coeff2=0.5079548825243172;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.0021646851456075737;
    coeff2=0.048228667223859864;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=-0.02961357901813876;
    coeff2=0.18883891859269755;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1=0.972022432714026;
    coeff2=0.04433350924092875;
    Bb8_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.0014416692554755693;
    coeff2=0.5079548825243172;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00010342360002606436;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.00025042089061726836;
    Bb8_3= coeff1*Bb8_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0142963038451998;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0000472074218614;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.15997208123052736;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=7.905233109508188;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1172284755624928;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.04415770962086628;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02125094412156715;
    Bb8_4= coeff1*Bb8_3+coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=6.015377936297571e-8;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=6.158327593023553e-7;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.05303925642383531;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0223309680940942;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.13164716678567034;
    Bb8_5= coeff1*Bb8_4+coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=2.9360008912888267e-7;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.01941558085813326;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.05424050332180066;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0995889901256144;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Ba5_2 with operation: lincomb
    coeff1=-0.0022569362567348168;
    coeff2=0.0005186821886494571;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.018857771843918;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0232717682963619;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.08242762761912202;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb5_2 with operation: lincomb
    coeff1=-0.016771544128867315;
    coeff2=0.0009494454604838552;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=7.889183859666258e-5;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0009930070661988107;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=-2.623470280367754e-5;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=4.3293107207049323e-7;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0061288549915095;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0486316092008641;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0031662589364607322;
    Bb8_6= coeff1*Bb8_5+coeff2*B5;
    % Computing Bb6_2 with operation: lincomb
    coeff1=0.020097870232187873;
    coeff2=0.05334311646938747;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.034593747722119206;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.07438835714613701;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.10144045604887632;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9872018647666524;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2841697393727225;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=4.056418623684063e-6;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.01581545457058901;
    Bb8_7= coeff1*Bb8_6+coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1=0.972022432714026;
    coeff2=0.04433350924092875;
    Ba8_2= coeff1*I+coeff2*A;
    % Computing Ba8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.00025042089061726836;
    Ba8_3= coeff1*Ba8_2+coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02125094412156715;
    Ba8_4= coeff1*Ba8_3+coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.13164716678567034;
    Ba8_5= coeff1*Ba8_4+coeff2*B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0031662589364607322;
    Ba8_6= coeff1*Ba8_5+coeff2*B5;
    % Computing Ba8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.01581545457058901;
    Ba8_7= coeff1*Ba8_6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=0.00883131198257824;
    coeff2=0.14508301662133044;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.07731260692174748;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.01972149511967523;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2595200313853995;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.9401871371954384;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9656190992569235;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9769076343982577;
    Ba8= coeff1*Ba8_7+coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9769076343982577;
    Bb8= coeff1*Bb8_7+coeff2*B7;
    % Computing B8 with operation: mult
    B8=Ba8 * Bb8;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=-5.923166377309273e-7;
    T2k8= coeff1*T2k7+coeff2*B7;
    % Computing T2k10 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0000001475289326;
    T2k10= coeff1*T2k8+coeff2*B8;
    output=T2k10;
end

