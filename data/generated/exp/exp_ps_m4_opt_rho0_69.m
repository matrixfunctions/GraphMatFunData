function output = exp_ps_m4_opt_rho0_69(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Bb4_3 Ba3_2 Ba3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 T2k5 Bb5 B5 T2k7
    % Computing Bb2 with operation: lincomb
    coeff1 = 1.7372701751396924;
    coeff2 = -0.21601444678956824;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = 0.2866026571168137;
    coeff2 = -0.06808693499430804;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = -0.6108235957535929;
    coeff2 = -0.014709720427431594;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = 1.7372701751396924;
    coeff2 = -0.21601444678956824;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.1770525298266241;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.020859019885585126;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 1.0508101696874834;
    coeff2 = 1.040160435318157;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.474210559924759;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -0.3880998628583675;
    coeff2 = 2.232643563510906;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.6627704342768925;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = -1.8244556096069469;
    coeff2 = 1.0277285724068361;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.2223238419486635;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1 = 0.37176853543253857;
    coeff2 = 5.363140388026199;
    T2k2 = coeff1*I + coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.0137337080611135;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 2.92870970078987;
    coeff2 = 1.395210200929235;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9208910780611814;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.745273367542917;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.002429240292310237;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.104740854605295;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5901519971324986;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.6019799776198416;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.7913945361070499;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -2.338782896508547;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6037393214424426;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.6773088417874666;
    T2k7 = coeff1*T2k5 + coeff2*B5;
    output = T2k7;
end

