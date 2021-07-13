function output=exp_ps_opt_rho2_7(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9
    % Computing Ba3_2 with operation: lincomb
    coeff1=-0.8558107709451741;
    coeff2=6.654420725242113;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba6_2 with operation: lincomb
    coeff1=0.0035871003621533335;
    coeff2=-0.199515898070094;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing T2k2 with operation: lincomb
    coeff1=-4.433844812389262;
    coeff2=0.14174985095473627;
    T2k2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=2.373110447113857;
    coeff2=-11.496160069358545;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.7086288697769633;
    coeff2=1.5171963077831512;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=-0.7929434423854941;
    coeff2=5.6454298259935385;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=2.373110447113857;
    coeff2=-11.496160069358545;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-8.805953709057599;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.1464046050697995;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=3.3717237028855647;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=14.35374095270583;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=16.186403987022313;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=4.521013785174327e-9;
    coeff2=-0.016780854076753563;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.12296702237640218;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=0.01201228900944982;
    coeff2=-1.1065314556032247e-12;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-2.0779221004950836e-15;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing Ba7_2 with operation: lincomb
    coeff1=-0.09860921274946392;
    coeff2=0.2639396617162851;
    Ba7_2= coeff1*I+coeff2*A;
    % Computing Ba7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-4.745046748609644;
    Ba7_3= coeff1*Ba7_2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=-0.13628986461873388;
    coeff2=-0.027077600537154935;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.004774676777784938;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=2.6186718329845604;
    coeff2=4.547282965965041;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=5.695287571427549;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=25.592199199005485;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.6370814945587178;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.012069483753455755;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.1554804561626417;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.914091240339738;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.009542640674472407;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=-16.423040024910904;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.6798059142303008;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-5.249610150102727e-19;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=4.267589684021346e-27;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=3.7664383775767242e-6;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8690669836434851;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.0770978480112308;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=2.1074592605382616;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.5282957524727887;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.04873878019200554;
    T2k7= coeff1*T2k6+coeff2*B6;
    % Computing Bb7_2 with operation: lincomb
    coeff1=5.692302453618037;
    coeff2=-2.4680601957893162;
    Bb7_2= coeff1*I+coeff2*A;
    % Computing Bb7_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7182801012293478;
    Bb7_3= coeff1*Bb7_2+coeff2*B2;
    % Computing Bb7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.5030745298888237;
    Bb7_4= coeff1*Bb7_3+coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1=1.0;
    coeff2=3.145388920252181;
    Ba7_4= coeff1*Ba7_3+coeff2*B3;
    % Computing Bb7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.2802466342872092;
    Bb7_5= coeff1*Bb7_4+coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1=1.0;
    coeff2=3.6815285908979076;
    Ba7_5= coeff1*Ba7_4+coeff2*B4;
    % Computing Bb7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-7.214590944431316e-5;
    Bb7_6= coeff1*Bb7_5+coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.00016028279238559438;
    Ba7_6= coeff1*Ba7_5+coeff2*B5;
    % Computing Ba7 with operation: lincomb
    coeff1=1.0;
    coeff2=-9.362127616405266e-5;
    Ba7= coeff1*Ba7_6+coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8753927048889629;
    Bb7= coeff1*Bb7_6+coeff2*B6;
    % Computing B7 with operation: mult
    B7=Ba7 * Bb7;
    % Computing T2k9 with operation: lincomb
    coeff1=1.0;
    coeff2=1.3559720150054728;
    T2k9= coeff1*T2k7+coeff2*B7;
    output=T2k9;
end

