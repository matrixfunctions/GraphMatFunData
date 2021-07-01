function output=denman_beavers_opt(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Bb4_3 Ba3_2 Ba3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 T2k5 Bb5 B5 T2k7
    % Computing Bb2 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.011141376535527443;
    Bb2= coeff1*I+coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5152503668289159;
    Ba4_2= coeff1*I+coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1=1.0;
    coeff2=0.03845623901890617;
    Bb4_2= coeff1*I+coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1=1.0;
    coeff2=0.7610413081657074;
    Ba2= coeff1*I+coeff2*A;
    % Computing B2 with operation: ldiv
    B2=Ba2 \ Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.10150239642416546;
    Ba4_3= coeff1*Ba4_2+coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.026744368550352622;
    Bb4_3= coeff1*Bb4_2+coeff2*B2;
    % Computing Ba3_2 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.07127247848266706;
    Ba3_2= coeff1*I+coeff2*A;
    % Computing Ba3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.9541397430184381;
    Ba3= coeff1*Ba3_2+coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1=0.49421282325538585;
    coeff2=0.05943417585870563;
    Ba5_2= coeff1*I+coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5017024742690791;
    Ba5_3= coeff1*Ba5_2+coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1=1.0;
    coeff2=0.12762852235945304;
    Bb5_2= coeff1*I+coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0015101375211018255;
    Bb5_3= coeff1*Bb5_2+coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1=0.28457285753903816;
    coeff2=0.11495283295141033;
    T2k2= coeff1*I+coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.0344930890633602;
    T2k3= coeff1*T2k2+coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.02575361014723645;
    Bb3_2= coeff1*I+coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.004282806676350282;
    Bb3= coeff1*Bb3_2+coeff2*B2;
    % Computing B3 with operation: ldiv
    B3=Ba3 \ Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.016950835467818528;
    Ba4= coeff1*Ba4_3+coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.029043605652931636;
    Bb4= coeff1*Bb4_3+coeff2*B3;
    % Computing B4 with operation: ldiv
    B4=Ba4 \ Bb4;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.5276435956563207;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    Ba5= coeff1*Ba5_4+coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0;
    coeff2=0.34606871085737384;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.14486471812413823;
    Bb5_4= coeff1*Bb5_3+coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0;
    coeff2=-0.006442438999085179;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0;
    coeff2=0.01583234295751802;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: ldiv
    B5=Ba5 \ Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1=1.0;
    coeff2=1.4670072325444539;
    T2k7= coeff1*T2k5+coeff2*B5;
    output=T2k7;
end

