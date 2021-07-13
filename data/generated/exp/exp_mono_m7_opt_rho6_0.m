function output=exp_mono_m7_opt_rho6_0(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10
    % Computing Ba3_2 with operation: lincomb
    coeff1=-0.03985250441225717;
    coeff2=0.1604494657928954;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=-0.00011974681141880265;
    coeff2=-9.609608204465964e-5;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=-0.00018687382272586462;
    coeff2=-0.00010155826657765157;
    T2k2= coeff1*I+coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1=-4.539313763978297e-6;
    coeff2=2.446796588580378e-6;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.3044412670842123;
    coeff2=0.4680830251883016;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.19853207576158147;
    coeff2=0.4568156227748025;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.01206538482409284;
    coeff2=-0.007524875465546422;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=-0.011223283513386698;
    coeff2=0.5004951105037614;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1=0.9834849370950302;
    coeff2=0.40790543324377204;
    Bb8_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.19853207576158147;
    coeff2=0.4568156227748025;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.029166175784732254;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.4202712650569836;
    Bb8_3= coeff1*Bb8_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.32310560219905493;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9304976822330905;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=4.215986081751763e-5;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0061270230013822065;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9639420128969182;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.004160420231960609;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.15990304681757225;
    Bb8_4= coeff1*Bb8_3+coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0001809250904584009;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.00011452178091543187;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0587301450366299;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0020530677449521447;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.033732264254459236;
    Bb8_5= coeff1*Bb8_4+coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=-3.544767352137871e-5;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.13025928619739066;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8471750409761895;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.06884543294512852;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Ba5_2 with operation: lincomb
    coeff1=-0.0014350882470920096;
    coeff2=-0.0029204003557249366;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.004322247579028231;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.013953282030764238;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8220608344548307;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.8268548885673269;
    coeff2=0.32692355515774696;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.05385438234818326;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.003348134917199928;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=4.728916919898747e-5;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.1868405154805229e-5;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8150363087806349;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.1203307382956932;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.008536371457893934;
    Bb8_6= coeff1*Bb8_5+coeff2*B5;
    % Computing Bb6_2 with operation: lincomb
    coeff1=-0.10687141220929966;
    coeff2=0.49573353466641706;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.09669233528231447;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.008584243736221898;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0002579652548561839;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=-6.616194274488254e-8;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9948104309724103;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.9506111159992416e-6;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0010322498247622353;
    Bb8_7= coeff1*Bb8_6+coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1=0.9834849370950302;
    coeff2=0.40790543324377204;
    Ba8_2= coeff1*I+coeff2*A;
    % Computing Ba8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.4202712650569836;
    Ba8_3= coeff1*Ba8_2+coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.15990304681757225;
    Ba8_4= coeff1*Ba8_3+coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.033732264254459236;
    Ba8_5= coeff1*Ba8_4+coeff2*B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.008536371457893934;
    Ba8_6= coeff1*Ba8_5+coeff2*B5;
    % Computing Ba8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0010322498247622353;
    Ba8_7= coeff1*Ba8_6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=-8.127184169897754e-8;
    coeff2=0.5009422778394859;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.10725741908100991;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.011203185515508135;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0005125194848333935;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.4968901291255363e-5;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.941362084939232e-8;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1=1.0;
    coeff2=4.413245234287325e-5;
    Ba8= coeff1*Ba8_7+coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1=1.0;
    coeff2=4.413245234287325e-5;
    Bb8= coeff1*Bb8_7+coeff2*B7;
    % Computing B8 with operation: mult
    B8=Ba8 * Bb8;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.0864858600469639e-7;
    T2k8= coeff1*T2k7+coeff2*B7;
    % Computing T2k10 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0004223592531323;
    T2k10= coeff1*T2k8+coeff2*B8;
    output=T2k10;
end

