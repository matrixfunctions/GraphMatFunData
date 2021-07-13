function output=exp_mono_m5_opt_rho1_68(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1=-4.6346414362462995e-5;
    coeff2=0.0015183205977868746;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.0010904030413732757;
    coeff2=1.0008760081747117;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-5.439593588786742e-5;
    coeff2=8.642508638357873e-5;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=0.0010337804923906489;
    coeff2=0.9970317188215052;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.0010904030413732757;
    coeff2=1.0008760081747117;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9996379335010811;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0011897338053029579;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.08229676673753394;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1=-1.2408057646229667e-6;
    coeff2=-3.1395014660350098e-6;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.8699534803824349e-6;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=-8.30025660582309e-6;
    coeff2=-1.614854194638885e-5;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.2536139794669188e-5;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.0007031814361206821;
    coeff2=1.0010882362282207;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.05182299699949662;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=0.9999994211561953;
    coeff2=0.9989210320213427;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.49825582022614506;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=-2.7188726824232142e-6;
    coeff2=0.9998021510172249;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.06065097803531471;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.001530127478825033;
    coeff2=0.9994108881899734;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.030492786462992225;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9997653650686295;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0031469743629754005;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.005091212395247093;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0010411131935665154;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9999723778369385;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00034385565143820764;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.16623132434744778;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.036597641936229344;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.007748443487686225;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=2.3223775362920453e-6;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0002736987094186782;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0007434150049329373;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.999999459268548;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=2.2549509122548797e-6;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.004190644602100282;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0004888481055895867;
    T2k8= coeff1*T2k6+coeff2*B6;
    output=T2k8;
end

