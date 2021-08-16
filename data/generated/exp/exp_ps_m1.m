function output = exp_ps_m1(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: B_0_1 A2 P0
    % Computing B_0_1 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    B_0_1 = coeff1*I + coeff2*A;
    % Computing A2 with operation: mult
    A2 = A * A;
    % Computing P0 with operation: lincomb
    coeff1 = 0.5;
    coeff2 = 1.0;
    P0 = coeff1*A2 + coeff2*B_0_1;
    output = P0;
end

