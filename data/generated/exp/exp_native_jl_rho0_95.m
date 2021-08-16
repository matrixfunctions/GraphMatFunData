function output = exp_native_jl_rho0_95(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: A2 Ua2 V2 A4 Ua3 V3 A6 Ua U V Z X P
    % Computing A2 with operation: mult
    A2 = A * A;
    % Computing Ua2 with operation: lincomb
    coeff1 = 8.64864e6;
    coeff2 = 277200.0;
    Ua2 = coeff1*I + coeff2*A2;
    % Computing V2 with operation: lincomb
    coeff1 = 1.729728e7;
    coeff2 = 1.99584e6;
    V2 = coeff1*I + coeff2*A2;
    % Computing A4 with operation: mult
    A4 = A2 * A2;
    % Computing Ua3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1512.0;
    Ua3 = coeff1*Ua2 + coeff2*A4;
    % Computing V3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 25200.0;
    V3 = coeff1*V2 + coeff2*A4;
    % Computing A6 with operation: mult
    A6 = A2 * A4;
    % Computing Ua with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Ua = coeff1*Ua3 + coeff2*A6;
    % Computing U with operation: mult
    U = Ua * A;
    % Computing V with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 56.0;
    V = coeff1*V3 + coeff2*A6;
    % Computing Z with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -1.0;
    Z = coeff1*V + coeff2*U;
    % Computing X with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    X = coeff1*V + coeff2*U;
    % Computing P with operation: ldiv
    P = Z \ X;
    output = P;
end

