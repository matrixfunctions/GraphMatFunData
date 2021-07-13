function output=exp_sastre_m3(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: T2k2 B2 T2k3 Bb3 Ba4_3 B3 Ba4 Bb4 B4 T2k4 T2k6
    % Computing T2k2 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb3 with operation: lincomb
    coeff1=0.01992047682223989;
    coeff2=0.004980119205559973;
    Bb3= coeff1*A+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=0.8765009801785554;
    coeff2=0.07665265321119147;
    Ba4_3= coeff1*A+coeff2*B2;
    % Computing B3 with operation: mult
    B3=B2 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=0.1225521150112075;
    coeff2=1.0;
    Bb4= coeff1*B2+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=2.974307204847627;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k6= coeff1*T2k4+coeff2*B4;
    output=T2k6;
end

