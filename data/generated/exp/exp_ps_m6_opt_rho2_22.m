function output=exp_ps_m6_opt_rho2_22(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1=0.14883975012114015;
    coeff2=0.8408521785457064;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=2.880470058346373e-5;
    coeff2=-0.04253750167360615;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=0.7426007763589999;
    coeff2=0.23958838183794046;
    T2k2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.6932267294704894;
    coeff2=0.42019701808349896;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.03430490914845939;
    coeff2=0.3596800085189574;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=1.0264634378476645;
    coeff2=0.6572335936947115;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.6932267294704894;
    coeff2=0.42019701808349896;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8532569213076548;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.6441812976652185;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.13303235916009215;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2456652487277263;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.054585888069760787;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=1.0499238030144577e-10;
    coeff2=0.0012477776392267876;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.000784769897476988;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=4.054924703391629e-7;
    coeff2=-2.3045432301617836e-8;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-7.77392066692997e-10;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1=0.9986923999012539;
    coeff2=0.04556161788104884;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.39768438765561237;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=0.0004535421655560562;
    coeff2=2.887116538841799e-5;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=2.7169875443715466e-5;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=1.042613839944681;
    coeff2=0.3549827751063252;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.19349580056392088;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9021858527038544;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.028297918922423518;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.1358540783808526e-6;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.10903345552815719;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9996647988764669;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.232855953231344e-7;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.19442822324865222;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.12981790073255753;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=3.2001904605874257e-10;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.1346120403005061e-12;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=2.3660170281050782e-7;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0000000017458077;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8152592432110274;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9444261197925446;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.7522893501377047e-7;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0013719348524986168;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=-0.03834000555185737;
    coeff2=0.014143811582239133;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.008028642819640242;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0023561023048899606;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.22180315416394053;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=3.575682812162475e-5;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.094529320861645;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=4.24521104018061e-5;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.9466631013553095e-7;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00024920094919997024;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0000015155621724;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9988773313640329;
    T2k9= coeff1*T2k7+coeff2*B7;
    output=T2k9;
end

