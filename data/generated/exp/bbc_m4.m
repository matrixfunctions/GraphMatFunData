function output=bbc_m4(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba5_2 Bb5_2 T2k2 B2 T2k3 Ba5_3 Bb5_3 Ba4_3 Bb4_3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 Bb5 B5 T2k7
    % Computing Ba5_2 with operation: lincomb
    coeff1=4.811693118299809;
    coeff2=1.1510994882542136;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.21169311829980944;
    coeff2=0.15822438471572672;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=-0.018602320514620553;
    coeff2=-0.005007023225733177;
    T2k2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.5734201229605222;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03318960838392776;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.1656351694367274;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=-0.13181061013830184;
    coeff2=-0.02027855540589259;
    Ba4_3= coeff1*A+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=-0.13181061013830184;
    coeff2=-0.02027855540589259;
    Bb4_3= coeff1*A+coeff2*B2;
    % Computing B3 with operation: mult
    B3=B2 * A;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.006759518468630863;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.006759518468630863;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.012516177931579242;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.13339969394389206;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.010786277931579243;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    T2k7= coeff1*T2k4+coeff2*B5;
    output=T2k7;
end

