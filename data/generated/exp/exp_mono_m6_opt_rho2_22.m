function output = exp_mono_m6_opt_rho2_22(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 0.01492095365794394;
    coeff2 = 0.4548827343944871;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = -9.936953530550332e-5;
    coeff2 = 0.0004949489851157203;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 0.9461427385539992;
    coeff2 = 0.6846781699479447;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = 0.3606683467582314;
    coeff2 = 0.9022563141368555;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -0.00010283625521797217;
    coeff2 = 0.04865972405992692;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 0.20434201343253855;
    coeff2 = 1.0771903272919514;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = 0.3606683467582314;
    coeff2 = 0.9022563141368555;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.40742189739028956;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9805199719779003;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.005458404492452591;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.3428623107685864;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1042620847957364;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -0.0002492898506134758;
    coeff2 = -0.002520543826889288;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.19205783802229012;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = -0.011182019398311033;
    coeff2 = 1.0828912368373977;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1123228135458941;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = 0.2015415372769708;
    coeff2 = -1.2056970025906528e-5;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0010961136609599802;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = -0.7988605479429666;
    coeff2 = 0.8780429466509988;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.12320927828592863;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 0.48747544345380134;
    coeff2 = 0.9573248464890546;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.07129556473752731;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0458479561631584;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0011131317329747153;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.006251856663919193;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.14428312267404939;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0868841114600982;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0001199235362706937;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.12217177615335534;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.022127837028090097;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0013519953360161272;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0003247638741731967;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0026742858471557607;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.047957137450033e-7;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1389733413094104;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.039861297957641054;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9899602421104243;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0002497865786489251;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = -2.4273343771322697e-7;
    coeff2 = 0.9970474322309476;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.09407776630184815;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.004953345778869348;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.007943900842528887;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00013153151535746322;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.12239957112110682;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.422001410109185e-7;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.7911209172110849;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9445221704805566;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.5631701775215771e-9;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 6.24719462909587e-6;
    T2k9 = coeff1*T2k7 + coeff2*B7;
    output = T2k9;
end

