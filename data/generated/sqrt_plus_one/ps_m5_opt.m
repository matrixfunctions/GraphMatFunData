function output=ps_m5_opt(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8
    % Computing Ba3_2 with operation: lincomb
    coeff1=-1.1763384175098646;
    coeff2=0.14819920360945626;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Bb2 with operation: lincomb
    coeff1=-0.15614734226906415;
    coeff2=1.0138950449120168;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=-0.32868543911603565;
    coeff2=0.04948779256327628;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=0.2715485387807492;
    coeff2=0.4515173495475029;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=-0.1561459961705319;
    coeff2=1.0138962847063822;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=Ba2 * Bb2;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7004279880959638;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=1.1645163087525163;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.45819098373171513;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba6_2 with operation: lincomb
    coeff1=0.7570005286093737;
    coeff2=0.007598464345057781;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.8669783274817463;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=0.1784350608747074;
    coeff2=-0.9109875948763699;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.20063142586799526;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=-0.9977699469177798;
    coeff2=0.7154077419527756;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.0589444472962166;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=0.2989395453295784;
    coeff2=2.080346787066744;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.5251262809722197;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb6_2 with operation: lincomb
    coeff1=0.6207474625776077;
    coeff2=-0.8881294075403886;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.24425285144993336;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=-0.15896658336808092;
    coeff2=1.03474409787952;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.6899181523788153;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: mult
    B3=Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.33379574828605685;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5179573824690393;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=Ba4 * Bb4;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.9312258295819543;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.6463200998393868;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8066435555717062;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-1.0814610150998485;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.004244973897574398;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0326540011916379;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.32574547503223705;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.14407160251709974;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.04488529740138668;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.6983909370796325;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.3011465317844473;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0407649022964955;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k6 with operation: lincomb
    coeff1=1.0;
    coeff2=0.8226732267898158;
    T2k6= coeff1*T2k5+coeff2*B5;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7171476120487754;
    T2k8= coeff1*T2k6+coeff2*B6;
    output=T2k8;
end

