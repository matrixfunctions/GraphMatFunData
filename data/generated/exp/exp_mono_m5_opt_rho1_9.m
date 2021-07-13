function output=exp_mono_m5_opt_rho1_9(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1=1.6001548380533195;
    coeff2=-0.9896362789001155;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=2.925822675582923;
    coeff2=0.2429580790486332;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=0.11062999124786087;
    coeff2=0.5847339177180232;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=3.2853771337963456;
    coeff2=0.3808592456926523;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=2.925822675582923;
    coeff2=0.2429580790486332;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=2.535123243045754;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.5630308051212967;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.5832562723145828;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1=-0.004633918773394507;
    coeff2=0.03077253772413558;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0048951603098945345;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=0.019638592876749326;
    coeff2=-0.06879708908316624;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.38255885470355305;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=1.7177215564670094;
    coeff2=1.356661682335682;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.283038970541202;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=1.089645716755357;
    coeff2=-0.05894345454785517;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.9565192611905264;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=0.00017310181177139793;
    coeff2=0.7091796380446608;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.8509852584224694;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=2.477932656299724;
    coeff2=0.12561938455295804;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.39676277541410987;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.4346195499890089;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.01314523058371281;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.029768444198531262;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.195010535882028;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.3539321663418642;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.002161503231210118;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.07941140929885598;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.06299059278775863;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.04857205018560095;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.39705603974476195;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0010901570728226485;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.6636659134668268;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.283915906769389;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=2.4473617044179496e-5;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.022010502754629162;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0020981712432435614;
    T2k8= coeff1*T2k6+coeff2*B6;
    output=T2k8;
end

