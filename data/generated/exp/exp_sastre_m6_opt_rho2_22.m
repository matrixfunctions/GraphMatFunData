function output=exp_sastre_m6_opt_rho2_22(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1=-0.005681896024400657;
    coeff2=0.023301056977275802;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=-0.07225676735259205;
    coeff2=0.40637944910599555;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=0.9658865558277218;
    coeff2=0.021151469113934267;
    T2k2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.002597101166550988;
    coeff2=0.9907226316794635;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.0371131949784736;
    coeff2=0.4122623848215059;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=-0.0993170088559773;
    coeff2=0.1075843791249758;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.002597101166550988;
    coeff2=0.9907226316794635;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00982117020658891;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.999729943278336;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.044422614636220305;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.017762277056814803;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03444785392998929;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=0.000852548961418871;
    coeff2=0.012977435119527533;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9913750471409082;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=-0.008832765141592277;
    coeff2=0.0022901357546005413;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0003677844696610227;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1=3.144403235434427;
    coeff2=0.5115991466872888;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.1652954236470456;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=-0.06348900975427792;
    coeff2=0.10321951688742986;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.04045126277415917;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=-0.011734927258548565;
    coeff2=0.000831101913142186;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=9.825887307216484e-5;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9691387395826522;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0058330089499217;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.028287777007324065;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.007013193540886666;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.21625812151699;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.13721732464684955;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0005931703373916548;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.09923354639251344;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.008680829603225283;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0006838601561276326;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.005656023271494485;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0702928736641044;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.09025415645822;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.1540082683372147;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0657896264887965;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0161502730596768;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=0.0055528225167482785;
    coeff2=0.3634011920588302;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.09318615672033989;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.06811866746595016;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=7.924897961202367;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.16962731509782766;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9723882198958138;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.85596241477658;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0692074219098335;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.009097568257407438;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0800725044535462;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9681060677762616;
    T2k9= coeff1*T2k7+coeff2*B7;
    output=T2k9;
end

