function output=exp_sastre_m7_opt_rho6_0(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10
    % Computing Ba3_2 with operation: lincomb
    coeff1=0.0015678308682144713;
    coeff2=-0.001644406122563135;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=0.012857678730185047;
    coeff2=0.20416569487077543;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=-0.004180474678744686;
    coeff2=-0.002570469257093357;
    T2k2= coeff1*I+coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1=3.1087497397974926;
    coeff2=0.3694214555352366;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1=-0.009556418198624978;
    coeff2=0.0008384641570218439;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.0015638303136970531;
    coeff2=0.5131399756353182;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.0025650924538428376;
    coeff2=0.0751072955748711;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=-0.029899374961962634;
    coeff2=0.17661587215939123;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1=0.9717659119346435;
    coeff2=0.03708072735488672;
    Bb8_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.0015638303136970531;
    coeff2=0.5131399756353182;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.00012764063512709018;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0032737750163317585;
    Bb8_3= coeff1*Bb8_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.013533154589413957;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0000648973491737;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.15935169296586307;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=7.9059413858030805;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1136335440836176;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.04205470017452324;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.021723450159739942;
    Bb8_4= coeff1*Bb8_3+coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0009454274909569563;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.003167848437529552;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.05543256972765084;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0216298885266226;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.13300038220488442;
    Bb8_5= coeff1*Bb8_4+coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.003321896755863846;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.01518879866885416;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.044060889132843686;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1026236327353052;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Ba5_2 with operation: lincomb
    coeff1=-0.0021189804037908145;
    coeff2=0.001037076938904245;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0216723193069674;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.022323328476928236;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.052946917852281214;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb5_2 with operation: lincomb
    coeff1=-0.016491731477218847;
    coeff2=0.0008905549732892586;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=9.855583452895887e-5;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0008621015249925553;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=-2.1080367849877218e-5;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.015009539610819524;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0071150699762386;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.048941742425579186;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.005257171591711241;
    Bb8_6= coeff1*Bb8_5+coeff2*B5;
    % Computing Bb6_2 with operation: lincomb
    coeff1=0.017396334722190973;
    coeff2=0.04276963803433957;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03219045127574837;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.07421310609041179;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.09089674435798574;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9929884164279147;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.24374173790191617;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.02988218764945254;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02345129553527665;
    Bb8_7= coeff1*Bb8_6+coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1=0.9717659119346435;
    coeff2=0.03708072735488672;
    Ba8_2= coeff1*I+coeff2*A;
    % Computing Ba8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0032737750163317585;
    Ba8_3= coeff1*Ba8_2+coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.021723450159739942;
    Ba8_4= coeff1*Ba8_3+coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.13300038220488442;
    Ba8_5= coeff1*Ba8_4+coeff2*B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.005257171591711241;
    Ba8_6= coeff1*Ba8_5+coeff2*B5;
    % Computing Ba8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02345129553527665;
    Ba8_7= coeff1*Ba8_6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=0.009133614148579657;
    coeff2=0.14669027898511353;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0775106388395176;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.020068924005133555;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2544089375564568;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.9371031112856034;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9819279402861139;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9810542134903524;
    Ba8= coeff1*Ba8_7+coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9810542134903524;
    Bb8= coeff1*Bb8_7+coeff2*B7;
    % Computing B8 with operation: mult
    B8=Ba8 * Bb8;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.004141999495986352;
    T2k8= coeff1*T2k7+coeff2*B7;
    % Computing T2k10 with operation: lincomb
    coeff1=1.0;
    coeff2=1.002949767210698;
    T2k10= coeff1*T2k8+coeff2*B8;
    output=T2k10;
end

