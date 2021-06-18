function output=sastre_m8(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba9_2 T2k2 Ba8_2 B2 Ba8_3 T2k3 Ba9_3 Bb6_3 Ba7_3 B3 Bb6_4 Ba9_4 Bb7_4 T2k4 Ba8_4 Ba7_4 B4 Ba8_5 Ba9_5 T2k5 Bb6_5 Bb7_5 Ba7_5 B5 Ba8_6 Bb7_6 Ba7_6 Bb6 B6 Ba8_7 Ba7 Bb7 B7 Ba8 B8 Ba9 B9 T2k11
    % Computing Ba9_2 with operation: lincomb
    coeff1=0.008333333333333333;
    coeff2=0.001388888888888889;
    Ba9_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k2= coeff1*I+coeff2*A;
    % Computing Ba8_2 with operation: lincomb
    coeff1=2.755731922398589e-7;
    coeff2=2.505210838544172e-8;
    Ba8_2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing Ba8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=2.08767569878681e-9;
    Ba8_3= coeff1*Ba8_2+coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Ba9_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0001984126984126984;
    Ba9_3= coeff1*Ba9_2+coeff2*B2;
    % Computing Bb6_3 with operation: lincomb
    coeff1=-1.023660713518307e-11;
    coeff2=-4.508311519886735e-13;
    Bb6_3= coeff1*A+coeff2*B2;
    % Computing Ba7_3 with operation: lincomb
    coeff1=-5.893435534477677e-5;
    coeff2=-3.013961104055248e-6;
    Ba7_3= coeff1*A+coeff2*B2;
    % Computing B3 with operation: mult
    B3=A * B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.980157255925737e-14;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba9_4 with operation: lincomb
    coeff1=1.0;
    coeff2=2.48015873015873e-5;
    Ba9_4= coeff1*Ba9_3+coeff2*B3;
    % Computing Bb7_4 with operation: lincomb
    coeff1=-5.100472475630675e-7;
    coeff2=-4.032817333361947e-8;
    Bb7_4= coeff1*B2+coeff2*B3;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.16666666666666666;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Ba8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.30531132637709e-10;
    Ba8_4= coeff1*Ba8_3+coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.502070379373464e-7;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=A * B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=7.55676813469492e-12;
    Ba8_5= coeff1*Ba8_4+coeff2*B4;
    % Computing Ba9_5 with operation: lincomb
    coeff1=1.0;
    coeff2=2.7557319223985893e-6;
    Ba9_5= coeff1*Ba9_4+coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.041666666666666664;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-9.210033748491798e-16;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-2.785084196756015e-9;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-6.770221628797445e-9;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=A * B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=4.024189993755686e-13;
    Ba8_6= coeff1*Ba8_5+coeff2*B5;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-3.294026127901678e-10;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.227011356117036e-10;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=-6.140022498994532e-17;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=B5 * Bb6;
    % Computing Ba8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.001023463999572971;
    Ba8_7= coeff1*Ba8_6+coeff2*B6;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba8= coeff1*Ba8_7+coeff2*B7;
    % Computing B8 with operation: mult
    B8=Ba8 * B5;
    % Computing Ba9 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba9= coeff1*Ba9_5+coeff2*B8;
    % Computing B9 with operation: mult
    B9=Ba9 * B5;
    % Computing T2k11 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k11= coeff1*T2k5+coeff2*B9;
    output=T2k11;
end

