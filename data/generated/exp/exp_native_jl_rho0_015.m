function output = exp_native_jl_rho0_015(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: A2 Ua U V Z X P
    % Computing A2 with operation: mult
    A2 = A * A;
    % Computing Ua with operation: lincomb
    coeff1 = 60.0;
    coeff2 = 1.0;
    Ua = coeff1*I + coeff2*A2;
    % Computing U with operation: mult
    U = Ua * A;
    % Computing V with operation: lincomb
    coeff1 = 120.0;
    coeff2 = 12.0;
    V = coeff1*I + coeff2*A2;
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

