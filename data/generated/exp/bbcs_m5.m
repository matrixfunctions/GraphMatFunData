function output=bbcs_m5(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: Ba6_2 Bb6_2 B2 Ba6_3 Bb6_3 T2k3 Ba5_3 B3 Bb6_4 Ba5_4 T2k4 Bb5_4 Ba6_4 B4 T2k5 Bb5 B5 Ba6_5 Bb6_5 Ba6 Bb6 B6 T2k8
    % Computing Ba6_2 with operation: lincomb
    coeff1=0.3420232802536553 + 0.0im;
    coeff2=0.0 + 1i*-0.2851997796332415;
    Ba6_2= coeff1*I+coeff2*A;
    % Computing Bb6_2 with operation: lincomb
    coeff1=2.9237775839655367 + 0.0im;
    coeff2=0.0 + 1i*1.4451330034748826;
    Bb6_2= coeff1*I+coeff2*A;
    % Computing B2 with operation: mult
    B2=A * A;
    % Computing Ba6_3 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=0.047347067331271094 + 0.0im;
    Ba6_3= coeff1*Ba6_2+coeff2*B2;
    % Computing Bb6_3 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=0.1240818356655045 + 0.0im;
    Bb6_3= coeff1*Bb6_2+coeff2*B2;
    % Computing T2k3 with operation: lincomb
    coeff1=0.0 + 1i*-0.6604084076077131;
    coeff2=-1.093022784715649 + 0.0im;
    T2k3= coeff1*A+coeff2*B2;
    % Computing Ba5_3 with operation: lincomb
    coeff1=0.12 + 0.0im;
    coeff2=-0.0 + 1i*-0.008774760968797039;
    Ba5_3= coeff1*A+coeff2*B2;
    % Computing B3 with operation: mult
    B3=B2 * A;
    % Computing Bb6_4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.0 + 1i*-0.019571570936427238;
    Bb6_4= coeff1*Bb6_3+coeff2*B3;
    % Computing Ba5_4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.0009784845352378095 + 0.0im;
    Ba5_4= coeff1*Ba5_3+coeff2*B3;
    % Computing T2k4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=0.0 + 1i*0.2537715581771087;
    T2k4= coeff1*T2k3+coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1=0.0 + 1i*-0.12395369585828313;
    coeff2=-0.011202694841085593 + 0.0im;
    Bb5_4= coeff1*B2+coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.0 + 1i*-0.022186600635366212;
    Ba6_4= coeff1*Ba6_3+coeff2*B3;
    % Computing B4 with operation: mult
    B4=B3 * B3;
    % Computing T2k5 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=0.0005437426743473122 + 0.0im;
    T2k5= coeff1*T2k4+coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-0.0 + 1i*-1.2367240538259895e-5;
    Bb5= coeff1*Bb5_4+coeff2*B4;
    % Computing B5 with operation: mult
    B5=Ba5_4 * Bb5;
    % Computing Ba6_5 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=-9.74758985615379e-6 + 0.0im;
    Ba6_5= coeff1*Ba6_4+coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=2.425253007433925e-5 + 0.0im;
    Bb6_5= coeff1*Bb6_4+coeff2*B4;
    % Computing Ba6 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=1.0 + 0.0im;
    Ba6= coeff1*Ba6_5+coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=1.0 + 0.0im;
    Bb6= coeff1*Bb6_5+coeff2*B5;
    % Computing B6 with operation: mult
    B6=Ba6 * Bb6;
    % Computing T2k8 with operation: lincomb
    coeff1=1.0 + 0.0im;
    coeff2=1.0 + 0.0im;
    T2k8= coeff1*T2k5+coeff2*B6;
    output=T2k8;
end

