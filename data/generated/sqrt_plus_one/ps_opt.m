function output=ps_opt(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Bb4_3 Ba3_2 Ba3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 T2k5 Bb5 B5 T2k7
    % Computing Bb2 with operation: lincomb
    coeff1=0.166574948798779;
    coeff2=0.8082854563020354;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=0.04715365554193093;
    coeff2=0.15418701337518767;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=0.3432126022454872;
    coeff2=-0.11246622173262046;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.166574948798779;
    coeff2=0.8082854563020354;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.3797925656125231;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.046141102900243;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba3_2 with operation: lincomb
    coeff1=-0.19049161291756245;
    coeff2=-0.2272898739592557;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0768814160832578;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=0.46601596472571843;
    coeff2=0.14688479817683717;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.008944212309158862;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=-0.3596687172426221;
    coeff2=0.18769186938950672;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.05912889403697017;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=1.07746933317159;
    coeff2=0.5003011141213493;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.10032929104108079;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=-0.08695140350807312;
    coeff2=0.8768559204410328;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.8635667598418555;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1521649188195786;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.21328425050346164;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9604216332394367;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.0249903299528245;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.2019777447551502;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2678212927140638;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5451094283480894;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1812063585204995;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5334069828310514;
    T2k7= coeff1*T2k5+coeff2*B5;
    output=T2k7;
end

