function output=sid_m8(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Bb9_2 Bb2 Ba9_2 Ba2 B2 Bb9_3 Ba9_3 Ba3 B3 Ba9_4 Bb9_4 Ba8_3 Ba8_4 Ba4 B4 Ba8_5 Ba9_5 Bb9_5 Ba5 B5 Bb9_6 Ba8_6 Ba9_6 Bb7_4 Bb7_5 Bb7_6 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba8_7 Bb7 Ba7_3 Bb8_3 Bb8_4 Ba7_4 Bb8_5 Ba7_5 Ba7_6 Bb8_6 Ba7 Bb8_7 B7 Ba8 B8 Ba9 Bb9 B9 T2k11
    % Computing Bb9_2 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5;
    Bb9_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.0;
    coeff2=0.5;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba9_2 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5;
    Ba9_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.0;
    coeff2=0.5;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Bb9_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5029302610017967;
    Bb9_3= coeff1*Bb9_2+coeff2*B2;
    % Computing Ba9_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5029302610017967;
    Ba9_3= coeff1*Ba9_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=0.5;
    coeff2=0.0;
    Ba3= coeff1*A+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * B2;
    % Computing Ba9_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.07705596948494946;
    Ba9_4= coeff1*Ba9_3+coeff2*B3;
    % Computing Bb9_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.07705596948494946;
    Bb9_4= coeff1*Bb9_3+coeff2*B3;
    % Computing Ba8_3 with operation: lincomb
    coeff1=5.55344699042941;
    coeff2=2.991654767354374;
    Ba8_3= coeff1*A+coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2857950268422422;
    Ba8_4= coeff1*Ba8_3+coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1=0.5;
    coeff2=0.0;
    Ba4= coeff1*A+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03005135891320298;
    Ba8_5= coeff1*Ba8_4+coeff2*B4;
    % Computing Ba9_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.004985549176118462;
    Ba9_5= coeff1*Ba9_4+coeff2*B4;
    % Computing Bb9_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.004985549176118462;
    Bb9_5= coeff1*Bb9_4+coeff2*B4;
    % Computing Ba5 with operation: lincomb
    coeff1=0.5;
    coeff2=0.0;
    Ba5= coeff1*A+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * B4;
    % Computing Bb9_6 with operation: lincomb
    coeff1=1.0;
    coeff2=6.263526066651383e-5;
    Bb9_6= coeff1*Bb9_5+coeff2*B5;
    % Computing Ba8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.002742336655922557;
    Ba8_6= coeff1*Ba8_5+coeff2*B5;
    % Computing Ba9_6 with operation: lincomb
    coeff1=1.0;
    coeff2=6.263526066651383e-5;
    Ba9_6= coeff1*Ba9_5+coeff2*B5;
    % Computing Bb7_4 with operation: lincomb
    coeff1=0.03306559506631931;
    coeff2=0.002630043177655382;
    Bb7_4= coeff1*B2+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0002100333647757715;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=3.7537107416419e-5;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Bb6_3 with operation: lincomb
    coeff1=6.568407108494315e-7;
    coeff2=6.204734935438909e-8;
    Bb6_3= coeff1*A+coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=2.957106114715868e-9;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.556371639324141e-10;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.556371639324141e-11;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=B5 * Bb6;
    % Computing Ba8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=61.75954247606858;
    Ba8_7= coeff1*Ba8_6+coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing Ba7_3 with operation: lincomb
    coeff1=2.9757926317530323;
    coeff2=0.4155284057336423;
    Ba7_3= coeff1*A+coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1=-0.000131911811116888;
    coeff2=0.008139086096860678;
    Bb8_3= coeff1*A+coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.001121744731945438;
    Bb8_4= coeff1*Bb8_3+coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02479095151834799;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Bb8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=9.027588625491207e-5;
    Bb8_5= coeff1*Bb8_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.001283057135586989;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=3.501669195497238e-5;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=8.572383602707347e-6;
    Bb8_6= coeff1*Bb8_5+coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Bb8_7= coeff1*Bb8_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba8= coeff1*Ba8_7+coeff2*B7;
    % Computing B8 with operation: mult
    B8=Ba8 * Bb8_7;
    % Computing Ba9 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba9= coeff1*Ba9_6+coeff2*B8;
    % Computing Bb9 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Bb9= coeff1*Bb9_6+coeff2*B8;
    % Computing B9 with operation: mult
    B9=Ba9 * Bb9;
    % Computing T2k11 with operation: lincomb
    coeff1=0.0;
    coeff2=1.0;
    T2k11= coeff1*B7+coeff2*B9;
    output=T2k11;
end

