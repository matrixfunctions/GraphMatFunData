function output=exp_bbc_m5_opt_rho1_68(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1=-0.5512395258998516;
    coeff2=1.0032658350841217;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=-0.5675095568692214;
    coeff2=1.0513589278272073;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.0001314758343281401;
    coeff2=-0.009487670689278912;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=-0.0001314758343281401;
    coeff2=-0.009487670689278912;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=-0.5675095568692214;
    coeff2=1.0513589278272073;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.016975592917509923;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.2972477049482875;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.2972477049482875;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1=-11.042102494511589;
    coeff2=1.9369906298602684;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.13170872238552217;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=0.03195710896762638;
    coeff2=-0.2728747582815884;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.014819300077579842;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=-0.02894640931162959;
    coeff2=0.06013102031320287;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.2658107563575296;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=-0.01092637958920494;
    coeff2=0.3494106318836958;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.4203994390309589;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=-0.03815007289738427;
    coeff2=-0.16616807134202866;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.049988314856347284;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=-0.0484315404850608;
    coeff2=-0.5225875384324369;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0136434771969927;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9648816429267677;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9648816429267677;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.04289915281255721;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0011708512130180722;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=-4.637252920584828e-8;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.393714904854506e-5;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8132020643417125;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0005123761117618341;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0019952701863870247;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.03188188105549395;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=-3.909390915989026e-6;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=2.020645293993015e-5;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8068418040933824;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1923960333756067;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.009154750940822942;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1955049607512138;
    T2k8= coeff1*T2k6+coeff2*B6;
    output=T2k8;
end

