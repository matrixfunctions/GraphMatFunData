function output=sastre_m6(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba7_2 B2 Ba7_3 Bb3 B3 Ba7_4 Bb4 Ba4_3 Ba4 B4 Ba7_5 Ba6_3 Bb7_3 Bb5_3 B5 Ba6 Bb7_6 Bb6 B6 Bb7 B7 T2k9
    % Computing Ba7_2 with operation: lincomb
    coeff1=3.09646797193604;
    coeff2=0.772360321294401;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.1673139636901279;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Bb3 with operation: lincomb
    coeff1=0.001052151783051235;
    coeff2=0.0004675683454147702;
    Bb3= coeff1*A+coeff2*B2;
    % Computing B3 with operation: mult
    B3=B2 * Bb3;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=7.922322450524197;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=0.05317514832355802;
    coeff2=1.0;
    Bb4= coeff1*B2+coeff2*B3;
    % Computing Ba4_3 with operation: lincomb
    coeff1=0.2868706220817633;
    coeff2=-0.03289442879547955;
    Ba4_3= coeff1*A+coeff2*B2;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Ba6_3 with operation: lincomb
    coeff1=0.39689859154115;
    coeff2=0.02219811707032801;
    Ba6_3= coeff1*A+coeff2*B2;
    % Computing Bb7_3 with operation: lincomb
    coeff1=0.3229486011362678;
    coeff2=0.08092036376147299;
    Bb7_3= coeff1*A+coeff2*B2;
    % Computing Bb5_3 with operation: lincomb
    coeff1=0.002688394980266927;
    coeff2=0.0004675683454147702;
    Bb5_3= coeff1*A+coeff2*B2;
    % Computing B5 with operation: mult
    B5=B2 * Bb5_3;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba6= coeff1*Ba6_3+coeff2*B5;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.930814505527068;
    Bb7_6= coeff1*Bb7_3+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=0.0277140002806296;
    coeff2=1.0;
    Bb6= coeff1*B2+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7_5 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k9= coeff1*I+coeff2*B7;
    output=T2k9;
end

