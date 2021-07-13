function output=exp_mono_m4_opt_rho0_69(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Bb4_3 Ba3_2 Ba3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 T2k5 Bb5 B5 T2k7
    % Computing Bb2 with operation: lincomb
    coeff1=0.29569733081129235;
    coeff2=0.5579019761039045;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.10603592380569533;
    coeff2=0.05824938295333092;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=1.1108810217330913;
    coeff2=0.25686278900572646;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.29569733081129235;
    coeff2=0.5579019761039045;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5674218115639013;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.04714911885726358;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba3_2 with operation: lincomb
    coeff1=-0.22527364173809233;
    coeff2=0.12283893416894799;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7635839995489598;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=0.022799236610704456;
    coeff2=-1.0359744996702613;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.6837326299054575;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.015675034646166878;
    coeff2=1.260909852925862;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.3148486342463852;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=1.0216444204245163;
    coeff2=0.5251211761135689;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7800244667714745;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.8313380028114623;
    coeff2=0.1346607765677363;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.06933288385478588;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1157677193291784;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.10219516402901102;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.4707893269262639;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8155394030489035;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.13237659081921324;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.3364746558547452;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.48088503410167177;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.034247913899630375;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0674717167283821;
    T2k7= coeff1*T2k5+coeff2*B5;
    output=T2k7;
end

