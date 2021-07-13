function output=exp_sid_m5_opt_rho1_9(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1=-0.005962271096374199;
    coeff2=1.0033541057981434;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.029919432641542624;
    coeff2=0.9880527944917995;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=0.00036724475007954906;
    coeff2=-0.0001169342144972771;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=4.349957387597537e-5;
    coeff2=5.175987732045784e-5;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.029919432641542624;
    coeff2=0.9880527944917995;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.006633362317401653;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.00038475487324536156;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=3.984582078997789e-6;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1=0.4773034777980871;
    coeff2=0.2854265436738157;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.07274584364243558;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=-0.1355171122912882;
    coeff2=1.1473328323777334;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.08869410225976589;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=-0.06832215555028076;
    coeff2=0.03801394759424553;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.009069020035763865;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=0.9036674993687099;
    coeff2=1.0508277676937086;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.23923611287127894;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=0.22400344464331032;
    coeff2=0.7340989585819703;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.11736106888760392;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.0060878652738918255;
    coeff2=-0.006111731042202681;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0033794398005356;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9999998887371472;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=8.841956594385397e-7;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.027574521502167104;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.002958749892143305;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.2342100785171737;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=11.171900635837863;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.11879624303727956;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2616000497884793;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0033980901284630686;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.050826030283530055;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.3786103558350833;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.4285420863823937;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.229482143612308;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.545214800725622;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=3.414801497466826;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5760929600449584;
    T2k8= coeff1*T2k6+coeff2*B6;
    output=T2k8;
end

