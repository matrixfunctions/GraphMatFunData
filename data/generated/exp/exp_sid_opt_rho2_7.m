function output=exp_sid_opt_rho2_7(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1=0.013779624734115105;
    coeff2=1.0123617774427427;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=-0.0014340288130826606;
    coeff2=3.0214170368224798;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=0.9996354364827498;
    coeff2=1.0075667227584877;
    T2k2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.045435112444573475;
    coeff2=1.0223246692300103;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=0.0024179040551112;
    coeff2=0.9995075315498182;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=0.00023252325808871646;
    coeff2=4.612679859250955e-5;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.045435112444573475;
    coeff2=1.0223246692300103;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0417319729524033;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02096699267528694;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.28010545658135216;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.181554796011574;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.010778337993824194;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=7.788069195751279e-6;
    coeff2=8.668814613800501e-5;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0005851748274979261;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=4.367862801332661e-5;
    coeff2=4.771005078802816e-6;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=3.057699543578747e-7;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1=0.0003158514159932758;
    coeff2=8.290611574892452;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.7421134180734874;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=0.0009787233688326564;
    coeff2=-0.096167679660017;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.010643200832947241;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.0025872950992219526;
    coeff2=0.013669635797719902;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0123048626491173;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0010158820545039565;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9999831553043256;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.00024824310516588245;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0003724118165732193;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.000002894879292;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00012514140616432543;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.017961596596304457;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0006512172198766862;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.4369102577037629e-8;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=2.8715019361588243e-10;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=6.761523475937711e-6;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.228607675289351;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.017678291021512645;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0006690728448219196;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7423612498678709;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.036931956487199616;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=-0.0026162370361472948;
    coeff2=0.030992248965416677;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.009171582602004843;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0014654665589331566;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.30262885192060646;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=8.173851230520504e-5;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02034868019711193;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.119040499126371;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=249.89631597753709;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.227713124647877;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.4271746834360811e-5;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1186848822520254;
    T2k9= coeff1*T2k7+coeff2*B7;
    output=T2k9;
end

