function output = exp_bbc_m5(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Ba6_2 Bb6_2 B2 Ba6_3 Bb6_3 T2k3 Ba5_3 B3 Bb6_4 Ba5_4 T2k4 Bb5_4 Ba6_4 B4 T2k5 Bb5 B5 Ba6_5 Bb6_5 Ba6 Bb6 B6 T2k8
    % Computing Ba6_2 with operation: lincomb
    coeff1 = -11.058071288535288;
    coeff2 = 1.6125176868819238;
    Ba6_2 = coeff1*I + coeff2*A;
    % Computing Bb6_2 with operation: lincomb
    coeff1 = -0.09043168323908106;
    coeff2 = -0.06764045190713819;
    Bb6_2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.12477411482493252;
    Ba6_3 = coeff1*Ba6_2 + coeff2*B2;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.06759613017704597;
    Bb6_3 = coeff1*Bb6_2 + coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1 = 0.3978497494996451;
    coeff2 = 1.3678377846041172;
    T2k3 = coeff1*A + coeff2*B2;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = -0.10036558103014462;
    coeff2 = -0.00802924648241157;
    Ba5_3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = A * B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.029555257042931552;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.00089213849804573;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.49828962252538267;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = -0.09233646193671186;
    coeff2 = -0.016936493900208172;
    Bb5_4 = coeff1*B2 + coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.02257315581805103;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = B3 * B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.0006378981945947233;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.400867981820361e-5;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5_4 * Bb5;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.9579475957000978e-5;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.391802575160607e-5;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing T2k8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k8 = coeff1*T2k5 + coeff2*B6;
    output = T2k8;
end

