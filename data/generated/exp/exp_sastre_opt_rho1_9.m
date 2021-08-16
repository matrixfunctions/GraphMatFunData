function output = exp_sastre_opt_rho1_9(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1 = 0.5431852019042509;
    coeff2 = 2.808044477000553;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1 = -4.479293476259415;
    coeff2 = -1.965591461532149;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = 3.513744518233991;
    coeff2 = 3.480716229698247;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 0.02417939999948854;
    coeff2 = 0.013482072269100465;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -4.479293476259415;
    coeff2 = -1.965591461532149;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.072468368540587;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -7.395064443323127;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0015522395764164978;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1 = -2.57802584497714;
    coeff2 = -1.9681270306651504;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.5172562598353871;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = -2.7514538554861674;
    coeff2 = 28.533085265354913;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -32.70207976001747;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 1.4649998087002614;
    coeff2 = -1.1631539168725962;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 18.10426755193074;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1 = 9.679532976243642;
    coeff2 = 3.9359833423253106;
    T2k2 = coeff1*I + coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 6.090302916482898;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = 2.7927443516334343;
    coeff2 = -2.5111253025765686e-6;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.1720013130685409e-7;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 1.8758020055089961;
    coeff2 = -4.077994868395093;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.5726777535092064;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.9135175595994687;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -2.310576048441552e-6;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.1593312631548876e-10;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 15.117957339910209;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -31.69337657641421;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.7906362162645103e-10;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -2.7031124428719058;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -2.542481119421314;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -13.187053051683513;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.2838459715474883;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 27.426929433347695;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0930303538906508;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.17367745193891393;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -5.055543522880249e-19;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.6620101377768384;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.4265640842763014;
    T2k8 = coeff1*T2k6 + coeff2*B6;
    output = T2k8;
end

