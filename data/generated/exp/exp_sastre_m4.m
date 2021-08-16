function output = exp_sastre_m4(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: T2k2 B2 T2k3 Ba5_3 Bb4_3 B3 Bb4 Ba5_4 T2k4 Bb5_4 B4 Ba5 T2k5 Bb5 B5 T2k7
    % Computing T2k2 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.3093238729699403;
    coeff2 = 0.1955094199013519;
    Ba5_3 = coeff1*A + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 0.002193172316532563;
    coeff2 = 0.0002741465395665704;
    Bb4_3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = A * B2;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.569108992776174e-5;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.01626158346315151;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1168293067115003;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 0.03806343180936604;
    coeff2 = 0.017732587443103232;
    Bb5_4 = coeff1*B2 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = B3 * Bb4;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.018851944498568;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k7 = coeff1*T2k5 + coeff2*B5;
    output = T2k7;
end

