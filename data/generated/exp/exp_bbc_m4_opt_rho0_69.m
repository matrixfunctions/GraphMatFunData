function output=exp_bbc_m4_opt_rho0_69(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Bb4_3 Ba3_2 Ba3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 T2k5 Bb5 B5 T2k7
    % Computing Bb2 with operation: lincomb
    coeff1=0.08830117968963386;
    coeff2=1.6185302637301089;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=0.07611513553783929;
    coeff2=-0.1833134793841606;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=0.07611513553783929;
    coeff2=-0.1833134793841606;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.08830117968963386;
    coeff2=1.6185302637301089;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.005277078018897659;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.005277078018897659;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba3_2 with operation: lincomb
    coeff1=-0.012459653416356822;
    coeff2=0.10918913720213257;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0031491738964604;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=5.151267490244333;
    coeff2=0.3931553522425758;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0022244347351211098;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.14452058311618007;
    coeff2=0.15667930242227487;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.20054402254365047;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=-0.11353832118393399;
    coeff2=0.06878048308436636;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.4102780317055799;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.1634190055195114;
    coeff2=0.9966322449217647;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02185426146644069;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0019621792320780053;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0019621792320780053;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.000393833045141096;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2966584239338796;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.19117716256815653;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.015735447774198674;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.39486010387414244;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.4175497249260038;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.4114856620223937;
    T2k7= coeff1*T2k5+coeff2*B5;
    output=T2k7;
end

