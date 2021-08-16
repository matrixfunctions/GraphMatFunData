function output = exp_sastre_m6_opt_rho2_7(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 0.00203407200644078;
    coeff2 = -0.004855370488778173;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = 0.015636244039123465;
    coeff2 = 0.36835598195824903;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1 = 0.9720224032451862;
    coeff2 = 0.08866505998167876;
    T2k2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = 0.0014414740602762641;
    coeff2 = 1.015909177999167;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -0.00216500404580565;
    coeff2 = 0.09646644408880453;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = -0.029613578842783986;
    coeff2 = 0.377669985575162;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = 0.0014414740602762641;
    coeff2 = 1.015909177999167;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.014294207683808999;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0000472141119134;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.04415666392494207;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.00025013219458514967;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.05304182885989784;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -0.0022569437385435;
    coeff2 = 0.0010373403739564698;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.018857479770905;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = -0.0167704319125257;
    coeff2 = 0.0018989172164601057;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 7.889567871549961e-5;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1 = 3.1043463304552805;
    coeff2 = 0.7528521412951402;
    Ba7_2 = coeff1*I + coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.159971928559395;
    Ba7_3 = coeff1*Ba7_2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 0.020097851537399642;
    coeff2 = 0.10667353947597173;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.034594555518270974;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = -0.0068121857554928915;
    coeff2 = 0.0018059085648684517;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00010342034415488765;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.1172277610690486;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.022330965215494;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.07438766844651382;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.02327178284493455;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.08243100849197117;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.10144129673946056;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.021251100930175442;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.13164782415883697;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0009931016008313642;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -2.623385669729595e-5;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0031660405732008064;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.987202584803334;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.019415313402916037;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.054241015390228235;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0061285566950167;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.015815360491620606;
    T2k7 = coeff1*T2k6 + coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1 = 0.008831439952196395;
    coeff2 = 0.29016677501097143;
    Bb7_2 = coeff1*I + coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.07731361344655992;
    Bb7_3 = coeff1*Bb7_2 + coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.01972134822932684;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 7.905233218278997;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.2595231983751941;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0995871761977756;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.9401869208092082;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.048631444054980266;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.2841742967744319;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9656182105520791;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9769074674369014;
    T2k9 = coeff1*T2k7 + coeff2*B7;
    output = T2k9;
end

