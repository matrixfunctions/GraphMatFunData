function output=exp_ps_m7_opt_rho6_0(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10
    % Computing Ba3_2 with operation: lincomb
    coeff1=0.05431906428723421;
    coeff2=0.9217267072847067;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=-1.1142498138003632e-9;
    coeff2=-0.00044073573263142074;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=0.7754984117964325;
    coeff2=0.3475417761570095;
    T2k2= coeff1*I+coeff2*A;
    % Computing Ba7_2 with operation: lincomb
    coeff1=0.00030937511267045;
    coeff2=-0.004992129336289725;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Bb3_2 with operation: lincomb
    coeff1=0.7792960542315325;
    coeff2=0.4594622982985042;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=0.40493955401846876;
    coeff2=0.4037226736739943;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.07101339405931287;
    coeff2=0.2183449794974525;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=0.7611160284006683;
    coeff2=0.14582635451246748;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Bb8_2 with operation: lincomb
    coeff1=0.1078899138880161;
    coeff2=0.03205343917230924;
    Bb8_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=0.40493955401846876;
    coeff2=0.4037226736739943;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.085965940313993;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.005656891024698354;
    Bb8_3= coeff1*Bb8_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.6325401240010657;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.253861883536394;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.10003129254607988;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.21265617105393061;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5077334065861822;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.009944137968525259;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Bb8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.010962105218802998;
    Bb8_4= coeff1*Bb8_3+coeff2*B3;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.4359396287990991;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.007610711188155261;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.01346186908673783;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.005570337889186752;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0013745451070322613;
    Bb8_5= coeff1*Bb8_4+coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.43170388767357853;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.05904148896574505;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.16671013087805137;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.2535622236995746;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Ba5_2 with operation: lincomb
    coeff1=0.01237456812732672;
    coeff2=-0.007417121502402069;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.36461901499271315;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7318998265697871;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.42311906620920164;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.5304750742911883;
    coeff2=0.5804911103105471;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.15347975441496023;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.02847544670127242;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.017703603503626326;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.037684845865425066;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9891460760326148;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9892639607059927;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.004668205384343422;
    Bb8_6= coeff1*Bb8_5+coeff2*B5;
    % Computing Bb6_2 with operation: lincomb
    coeff1=5.222157680049785e-8;
    coeff2=-2.3658628414902608e-9;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-2.6488325501827903e-10;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=6.060309152649078e-11;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0661525964254732e-12;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=8.967326795327424e-15;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0457211702940734e-7;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.5754616856674562e-7;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.152223044748429e-5;
    Bb8_7= coeff1*Bb8_6+coeff2*B6;
    % Computing Ba8_2 with operation: lincomb
    coeff1=0.990632734847538;
    coeff2=-0.037157343305066634;
    Ba8_2= coeff1*I+coeff2*A;
    % Computing Ba8_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.09692046328345978;
    Ba8_3= coeff1*Ba8_2+coeff2*B2;
    % Computing Ba8_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.4186089963980845;
    Ba8_4= coeff1*Ba8_3+coeff2*B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.6631289790906699;
    Ba8_5= coeff1*Ba8_4+coeff2*B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7614234275617816;
    Ba8_6= coeff1*Ba8_5+coeff2*B5;
    % Computing Ba8_7 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.005632757710081e-8;
    Ba8_7= coeff1*Ba8_6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=0.006312989983639836;
    coeff2=3.05349059290055e-5;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.2830919585643768e-5;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=3.7913704630241052e-6;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=2.3225691022612054e-6;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=9.290484746562112e-8;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0000000001559715;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0713903017943822;
    Ba8= coeff1*Ba8_7+coeff2*B7;
    % Computing Bb8 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0009187967337483;
    Bb8= coeff1*Bb8_7+coeff2*B7;
    % Computing B8 with operation: mult
    B8=Ba8 * Bb8;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.0007760687742691415;
    T2k8= coeff1*T2k7+coeff2*B7;
    % Computing T2k10 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0058204367206869;
    T2k10= coeff1*T2k8+coeff2*B8;
    output=T2k10;
end

