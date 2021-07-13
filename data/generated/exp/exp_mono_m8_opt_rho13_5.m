function output=exp_mono_m8_opt_rho13_5(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Ba4_2 Bb9_2 Bb4_2 Ba9_2 Bb8_2 Ba5_2 Bb5_2 Bb6_2 Ba8_2 Bb7_2 B2 Bb3 Bb7_3 Bb6_3 Bb8_3 Ba4_3 Ba3 B3 Bb7_4 Ba5_3 Ba5_4 Bb5_3 Ba7_3 Ba9_3 Ba7_4 Ba8_3 Ba4 Bb6_4 Ba6_3 Ba9_4 T2k3 Bb4_3 Bb9_3 Bb4 B4 T2k4 Bb7_5 Bb9_4 Bb9_5 Ba5 Ba9_5 T2k5 Bb6_5 Ba8_4 Ba8_5 Bb5_4 Bb5 B5 T2k6 Bb9_6 Ba8_6 Ba9_6 Bb7_6 Bb6 Ba6_4 Bb8_4 Bb8_5 Ba6_5 Ba7_5 Ba6 B6 Bb9_7 Ba8_7 Ba9_7 Bb7 T2k7 Ba7_6 Bb8_6 Ba7 Bb8_7 B7 Ba9_8 Bb9_8 Ba8 Bb8 B8 Ba9 T2k9 Bb9 B9 T2k11
    % Computing Ba3_2 with operation: lincomb
    coeff1=-0.003584018414055765;
    coeff2=0.16376260367565243;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=4.437908455542293e-5;
    coeff2=0.0003258145596094315;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=0.008776543046053789;
    coeff2=-1.3078408314728238e-5;
    T2k2= coeff1*I+coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1=-2.8201037195155804e-6;
    coeff2=4.8728643050646066e-6;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.30528232877254763;
    coeff2=0.23746768414771152;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=0.0005515398938947843;
    coeff2=0.02533967378579724;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb9_2 with operation: lincomb
    coeff1=-0.06540541282818779;
    coeff2=6.846710342127941e-5;
    Bb9_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=-0.008872840493334803;
    coeff2=0.25571769129339744;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba9_2 with operation: lincomb
    coeff1=0.1346603332119568;
    coeff2=-5.9198518169123926e-5;
    Ba9_2= coeff1*I+coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1=1.0002604649982476;
    coeff2=0.2422828857323093;
    Bb8_2= coeff1*I+coeff2*A;
    % Computing Ba5_2 with operation: lincomb
    coeff1=-0.0012630062267903097;
    coeff2=-0.0010962696861478022;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.8293727566445452;
    coeff2=0.1675033881979259;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb6_2 with operation: lincomb
    coeff1=-0.10275230374603168;
    coeff2=0.25648676493721656;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Ba8_2 with operation: lincomb
    coeff1=1.0002604649982476;
    coeff2=0.2422828857323093;
    Ba8_2= coeff1*I+coeff2*A;
    % Computing Bb7_2 with operation: lincomb
    coeff1=0.004650329563627054;
    coeff2=0.2598813154789629;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0023273769915510135;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.006290014251994485;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.005628536962635935;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.022001533431922154;
    Bb8_3= coeff1*Bb8_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.017067262331527198;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.04917059370539749;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.010576806709956765;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00022839563887310042;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.013673005502563394;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0030839993650848517;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=2.2723646529748017e-6;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Ba9_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-8.660244928553117e-6;
    Ba9_3= coeff1*Ba9_2+coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0061294802470251075;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Ba8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.022001533431922154;
    Ba8_3= coeff1*Ba8_2+coeff2*B2;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9638185168573307;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.008329230127484285;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00022081927999800417;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Ba9_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-8.731780915688725e-5;
    Ba9_4= coeff1*Ba9_3+coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-9.364829337290504e-7;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0037256194747560907;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Bb9_3 with operation: lincomb
    coeff1=1.0;
    coeff2=2.738382227960167e-6;
    Bb9_3= coeff1*Bb9_2+coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.002052240135641998;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-4.881187067247219e-6;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0003949349245182783;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Bb9_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-6.3730674335546956e-6;
    Bb9_4= coeff1*Bb9_3+coeff2*B3;
    % Computing Bb9_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-5.739005263969626e-6;
    Bb9_5= coeff1*Bb9_4+coeff2*B4;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8216219421030205;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Ba9_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-2.1711606510903095e-5;
    Ba9_5= coeff1*Ba9_4+coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=-6.561388347938557e-7;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00021583757875408033;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Ba8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.15933029209747282;
    Ba8_4= coeff1*Ba8_3+coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03184541272812379;
    Ba8_5= coeff1*Ba8_4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0031897875289976277;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=3.234759791150232e-5;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.0279110807751132e-7;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Bb9_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-2.157568588204766e-6;
    Bb9_6= coeff1*Bb9_5+coeff2*B5;
    % Computing Ba8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.007635550247315368;
    Ba8_6= coeff1*Ba8_5+coeff2*B5;
    % Computing Ba9_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-5.9723490061415135e-6;
    Ba9_6= coeff1*Ba9_5+coeff2*B5;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=8.856309878298438e-6;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.0968604718569746e-7;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.13039438638797307;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Bb8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.15933029209747282;
    Bb8_4= coeff1*Bb8_3+coeff2*B3;
    % Computing Bb8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03184541272812379;
    Bb8_5= coeff1*Bb8_4+coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8475672106088131;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0688513777099278;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8145974867536132;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing Bb9_7 with operation: lincomb
    coeff1=1.0;
    coeff2=-3.3636703984819904e-7;
    Bb9_7= coeff1*Bb9_6+coeff2*B6;
    % Computing Ba8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0008615033049494069;
    Ba8_7= coeff1*Ba8_6+coeff2*B6;
    % Computing Ba9_7 with operation: lincomb
    coeff1=1.0;
    coeff2=-7.70252692547925e-7;
    Ba9_7= coeff1*Ba9_6+coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=7.553245375703866e-9;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=-5.458270135394687e-9;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.12038470670165373;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.007635550247315368;
    Bb8_6= coeff1*Bb8_5+coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9948034699390058;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0008615033049494069;
    Bb8_7= coeff1*Bb8_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing Ba9_8 with operation: lincomb
    coeff1=1.0;
    coeff2=-3.1855167508881976e-8;
    Ba9_8= coeff1*Ba9_7+coeff2*B7;
    % Computing Bb9_8 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.5592924041261544e-8;
    Bb9_8= coeff1*Bb9_7+coeff2*B7;
    % Computing Ba8 with operation: lincomb
    coeff1=1.0;
    coeff2=3.134369854398072e-5;
    Ba8= coeff1*Ba8_7+coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1=1.0;
    coeff2=3.134369854398072e-5;
    Bb8= coeff1*Bb8_7+coeff2*B7;
    % Computing B8 with operation: mult
    B8=Ba8 * Bb8;
    % Computing Ba9 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0081137407871135;
    Ba9= coeff1*Ba9_8+coeff2*B8;
    % Computing T2k9 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.06726697857924235;
    T2k9= coeff1*T2k7+coeff2*B8;
    % Computing Bb9 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9946441317159572;
    Bb9= coeff1*Bb9_8+coeff2*B8;
    % Computing B9 with operation: mult
    B9=Ba9 * Bb9;
    % Computing T2k11 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9964850858948753;
    T2k11= coeff1*T2k9+coeff2*B9;
    output=T2k11;
end

