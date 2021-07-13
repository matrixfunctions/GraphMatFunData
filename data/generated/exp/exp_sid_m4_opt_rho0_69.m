function output=exp_sid_m4_opt_rho0_69(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Bb4_3 Ba3_2 Ba3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 T2k5 Bb5 B5 T2k7
    % Computing Bb2 with operation: lincomb
    coeff1=9.94823321550182e-5;
    coeff2=1.0000210350390846;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=4.022059101053733e-5;
    coeff2=0.4019191862379914;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=2.1672609863472577e-5;
    coeff2=6.626362532960593e-5;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=9.94823321550182e-5;
    coeff2=1.0000210350390846;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.00873832228961145;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03233860555846545;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba3_2 with operation: lincomb
    coeff1=8.745327870027237e-6;
    coeff2=-4.197092420395855e-8;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9999999958392435;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=0.00011854240934339541;
    coeff2=2.2246677805188155;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2615074307076249;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=-7.86695052350628e-5;
    coeff2=-0.04148569168401734;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.023305174055381638;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=0.9999999940884431;
    coeff2=0.9999715234150552;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5918961310449276;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=-1.0660766661003063e-5;
    coeff2=0.0029440745098486642;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0004021729193984898;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0000192435512933;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0000837252012718;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=5.768958229431865;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9990123769097652;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-3.2734161536753197;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.023059406187239697;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=10.408133461599245;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9998635489305508;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9998619329256802;
    T2k7= coeff1*T2k5+coeff2*B5;
    output=T2k7;
end

