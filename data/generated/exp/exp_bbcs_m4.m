function output=exp_bbcs_m4(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba5_2 Bb5_2 T2k2 B2 T2k3 Ba5_3 Bb5_3 Ba4_3 Bb4_3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 Bb5 B5 T2k7
    % Computing Ba5_2 with operation: lincomb
    coeff1=2.6958430691533257 + 0.0im;
    coeff2=0.0 + 1i*0.05272871327381115;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Bb5_2 with operation: lincomb
    coeff1=2.6958430691533257 + 0.0im;
    coeff2=-0.0 + 1i*-1.3591092616886926;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=-6.267569853502023 + 0.0im;
    coeff2=0.0 + 1i*2.521796947120981;
    T2k2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=0.05786296656487002 + 0.0im;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.09896214548845832 + 0.0im;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.09896214548845832 + 0.0im;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=0.0 + 1i*0.13340427306445612;
    coeff2=0.020226020298183107 + 0.0im;
    Ba4_3= coeff1*A+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=0.0 + 1i*0.13340427306445612;
    coeff2=0.020226020298183107 + 0.0im;
    Bb4_3= coeff1*A+coeff2*B2;
    % Computing B3 with operation: mult
    B3=B2 * A;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.0 + 1i*-0.00674638241111651;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.0 + 1i*-0.00674638241111651;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=0.0 + 1i*0.007295441446830946;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=1.0 + 0.0im;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.0 + 1i*-0.0776668640807187;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=0.0 + 1i*0.015964794632994668;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=1.0 + 0.0im;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=1.0 + 0.0im;
    T2k7= coeff1*T2k4+coeff2*B5;
    output=T2k7;
end

