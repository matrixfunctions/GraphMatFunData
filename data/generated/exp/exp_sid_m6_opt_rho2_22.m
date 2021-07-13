function output=exp_sid_m6_opt_rho2_22(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1=0.0027676995071721643;
    coeff2=0.9956823779291472;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=0.00024723298310785875;
    coeff2=2.9161339295593187;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=0.9999550769133566;
    coeff2=0.9942510551400959;
    T2k2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.013941257411227636;
    coeff2=1.031053526686218;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=0.00010080993148845571;
    coeff2=1.000386931651619;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=0.00013960932977845965;
    coeff2=0.00027073100561329446;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.013941257411227636;
    coeff2=1.031053526686218;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.014565003723968966;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.018699401264499645;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.19294236215103988;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.23154125334090794;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00011767015392433715;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=1.0034850244473076e-5;
    coeff2=6.234107699654157e-5;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.00010709695216530405;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=-3.420829259855909e-5;
    coeff2=1.2724236347474151e-5;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=7.614113854507912e-7;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1=9.183179455916261e-5;
    coeff2=8.291209653793697;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.7562193397284078;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=0.0061957449875606665;
    coeff2=0.022600844232367568;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.035261063716355324;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.0014050828079002;
    coeff2=0.0025755913013428494;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9956033412300013;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0005245897678761895;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9999937753420692;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00401546253180321;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0004429662682755111;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9999999359403395;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0005271750223542335;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.011469528624315358;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0005156243401819641;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=4.506303286916229e-8;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=3.3732271485767847e-9;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.00015881528950764778;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0628031901019233;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.012746128809542087;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.001070557957948411;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.91478561889395;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.022922533875140892;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=-0.00015133350603989786;
    coeff2=0.028231887322807712;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.012411739707662116;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0004764037694311965;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.17985908646581833;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00018495252800435267;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.020190024440301577;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.05795634102314;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=249.89674278520076;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0634238665988902;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.000721682415651672;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0576256651496891;
    T2k9= coeff1*T2k7+coeff2*B7;
    output=T2k9;
end

