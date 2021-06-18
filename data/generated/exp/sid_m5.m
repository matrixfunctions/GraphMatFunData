function output=sid_m5(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: T2k2 B2 T2k3 Ba6_3 Ba5_3 Bb4_3 Bb6_3 B3 Bb4 Bb6_4 Ba5_4 T2k4 Bb5_4 Ba6_4 B4 Ba5 T2k5 Bb5 B5 Ba6_5 Bb6_5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing T2k2 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2863243726334417;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1=0.3112216227982407;
    coeff2=0.0852839259083158;
    Ba6_3= coeff1*A+coeff2*B2;
    % Computing Ba5_3 with operation: lincomb
    coeff1=0.9418613214806352;
    coeff2=0.06974348269544424;
    Ba5_3= coeff1*A+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=5.374708803114821e-5;
    coeff2=4.50085273957301e-6;
    Bb4_3= coeff1*A+coeff2*B2;
    % Computing Bb6_3 with operation: lincomb
    coeff1=0.6865706355662834;
    coeff2=0.1392249143769798;
    Bb6_3= coeff1*A+coeff2*B2;
    % Computing B3 with operation: mult
    B3=A * B2;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.16165883444488e-6;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03151382711608315;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.002005403977292901;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.08776036732867759;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1=-0.007544837153586671;
    coeff2=0.002852960512714315;
    Bb5_4= coeff1*B2+coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0292447258748138;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=B3 * Bb4;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.18995526739487723;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.829773504500424;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=11.173624766438472;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=3.23337016308538;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k8= coeff1*T2k6+coeff2*B6;
    output=T2k8;
end

