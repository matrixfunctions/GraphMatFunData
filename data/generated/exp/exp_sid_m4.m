function output = exp_sid_m4(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: T2k2 B2 T2k3 Bb3 B3 T2k4 Bb4 Ba4_3 Ba4 B4 T2k5 Ba5_3 Bb5_3 Ba5_4 Bb5_4 Ba5 Bb5 B5 T2k7
    % Computing T2k2 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5918659857804601;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb3 with operation: lincomb
    coeff1 = 0.002945531440279683;
    coeff2 = 0.0004018761610201036;
    Bb3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = B2 * Bb3;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -3.2733920099600837;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 0.03230762888122312;
    coeff2 = 1.0;
    Bb4 = coeff1*B2 + coeff2*B3;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 0.4017568440673568;
    coeff2 = -0.008709066576837676;
    Ba4_3 = coeff1*A + coeff2*B2;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 10.408017352313541;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 2.224209172496374;
    coeff2 = 0.2614927977298117;
    Ba5_3 = coeff1*A + coeff2*B2;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = -0.04130276365929783;
    coeff2 = 0.02338576034271299;
    Bb5_3 = coeff1*A + coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.768988513026145;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.023373194047115575;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
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

