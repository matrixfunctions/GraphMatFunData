function output = exp_sastre_m5_opt_rho1_9_low_relnorm(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1 = -0.357355571021064;
    coeff2 = 0.4105933675776958;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = -0.31551655352795516;
    coeff2 = 0.4801289831900131;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = -0.004620908141128493;
    coeff2 = 0.055057556018866645;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 0.011326704451019963;
    coeff2 = 0.002770498744594938;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -0.31551655352795516;
    coeff2 = 0.4801289831900131;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.09852026651116912;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.351290300034291;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00045550004840354903;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = 1.11114318963272;
    coeff2 = 0.4876551310270585;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.3710914185617396;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = 0.11176905543585865;
    coeff2 = 0.5925762886270637;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.4621864853147797;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 0.1215251552594178;
    coeff2 = 0.2145595643740629;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.04592597180399795;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1 = -0.0818796601588853;
    coeff2 = 0.003683007236493971;
    T2k2 = coeff1*I + coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.008733812199626457;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 0.9070495211990466;
    coeff2 = 0.5443460967617421;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5323384854781406;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 0.053558315750694814;
    coeff2 = -0.22014016998509814;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.922091412659715;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9607057168526784;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.00012416495659100344;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.12934538808055288;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.038392790243209074;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.24031719666304735;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.16322755510199;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.005536448121469381;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.003015592881979151;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.021375614386974554;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.08701143872628835;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.512024337751819;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.180154142825969;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6495647957423022;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6107117158792567;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.01836716605825683;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9780133820655359;
    T2k8 = coeff1*T2k6 + coeff2*B6;
    output = T2k8;
end

