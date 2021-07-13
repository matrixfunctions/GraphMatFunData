function output=exp_ps_m7_opt_rho3_59(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10
    % Computing Ba3_2 with operation: lincomb
    coeff1=0.07870552092616893;
    coeff2=0.9039871669928049;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=-9.869136468294225e-10;
    coeff2=-7.512320521145497e-5;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=0.7726601458443281;
    coeff2=0.27826200797600853;
    T2k2= coeff1*I+coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1=0.00041041842903290026;
    coeff2=0.0029599104645556024;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.8599464743752632;
    coeff2=0.29073009539473676;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.40288943046645087;
    coeff2=0.3243718504289327;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.07403200916806128;
    coeff2=0.22059858577526156;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=0.7557993691382382;
    coeff2=0.16186458120340602;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1=0.1230391215271295;
    coeff2=0.038029956154608906;
    Bb8_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.40288943046645087;
    coeff2=0.3243718504289327;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0877838768781113;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.010801521452913285;
    Bb8_3= coeff1*Bb8_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.6031002164926261;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.3209056390860863;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03421823700975471;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.07655965971121473;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5359493824375996;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0044353505165274365;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.019497085762494887;
    Bb8_4= coeff1*Bb8_3+coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.3803045282099576;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.11581791434832714;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.024014443649892782;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.010967880705307256;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0036174776941494966;
    Bb8_5= coeff1*Bb8_4+coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.3407182508161634;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.020309346530828727;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.195679851580968;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.45442773831109645;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Ba5_2 with operation: lincomb
    coeff1=-0.00013699934593016893;
    coeff2=0.14529476087609983;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.3261057160948171;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.710963649841826;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.34521134017746485;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.5028563878759467;
    coeff2=0.5066127059710326;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.19548297448489585;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.04481784891151489;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.023615787302546792;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.07858151588080117;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9850416366804396;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9367257958392154;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0021685857628971537;
    Bb8_6= coeff1*Bb8_5+coeff2*B5;
    % Computing Bb6_2 with operation: lincomb
    coeff1=9.211440361111637e-7;
    coeff2=-1.3293424291101191e-8;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-3.064863890686237e-9;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=7.122641395116041e-10;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=7.64822032891094e-12;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.3796783760471413e-13;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=-4.0465312386473035e-8;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.6667149787880624e-7;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.2133490010511278e-5;
    Bb8_7= coeff1*Bb8_6+coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1=0.7868458282010143;
    coeff2=-0.044300593493478396;
    Ba8_2= coeff1*I+coeff2*A;
    % Computing Ba8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.09816190377505705;
    Ba8_3= coeff1*Ba8_2+coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.45590771888971265;
    Ba8_4= coeff1*Ba8_3+coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.6380495620724085;
    Ba8_5= coeff1*Ba8_4+coeff2*B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7701122107742214;
    Ba8_6= coeff1*Ba8_5+coeff2*B5;
    % Computing Ba8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=3.0616077790644636e-8;
    Ba8_7= coeff1*Ba8_6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=0.00660875329127505;
    coeff2=9.739631514054463e-5;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=5.134838062997168e-5;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.9444505682154752e-5;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=9.194990346464199e-6;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-8.932667565096493e-8;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0000000001479474;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.040362632774053595;
    Ba8= coeff1*Ba8_7+coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1=1.0;
    coeff2=1.002439604236121;
    Bb8= coeff1*Bb8_7+coeff2*B7;
    % Computing B8 with operation: mult
    B8=Ba8 * Bb8;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00042398516888089884;
    T2k8= coeff1*T2k7+coeff2*B7;
    % Computing T2k10 with operation: lincomb
    coeff1=1.0;
    coeff2=1.009492506685841;
    T2k10= coeff1*T2k8+coeff2*B8;
    output=T2k10;
end

