function output = exp_sid_m6(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: T2k2 B2 T2k3 Ba6_3 Bb7_3 Bb5_3 Ba7_3 B3 Bb6_4 Bb7_4 T2k4 Bb5_4 Ba6_4 Ba7_4 B4 T2k5 Bb5 Ba6_5 Bb6_5 Bb7_5 Ba7_5 B5 Ba6 Bb7_6 Ba7_6 Bb6 B6 Ba7 B7 T2k9
    % Computing T2k2 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1969779342112314;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Ba6_3 with operation: lincomb
    coeff1 = 2.865001388641538;
    coeff2 = 0.1952545843107103;
    Ba6_3 = coeff1*A + coeff2*B2;
    % Computing Bb7_3 with operation: lincomb
    coeff1 = 0.03655234395347475;
    coeff2 = 0.01606091400855144;
    Bb7_3 = coeff1*A + coeff2*B2;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 2.294895435403922e-5;
    coeff2 = 1.406952242413849e-6;
    Bb5_3 = coeff1*A + coeff2*B2;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 8.29008575139441;
    coeff2 = 1.739158441630994;
    Ba7_3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = A * B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 0.02721930992200371;
    coeff2 = 0.002547056607231984;
    Bb6_4 = coeff1*B2 + coeff2*B3;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0001758035313846159;
    Bb7_4 = coeff1*Bb7_3 + coeff2*B3;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.03005000525808178;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 9.379681616092325e-8;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing Ba6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.01430688980356062;
    Ba6_4 = coeff1*Ba6_3 + coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.1965098904519709;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = A * B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.002243394407902074;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.172460202011541e-8;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing Ba6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.002024281516007681;
    Ba6_5 = coeff1*Ba6_4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.001204349003694297;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0002919349464582001;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.02018492049443954;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = B4 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Ba6 = coeff1*Ba6_5 + coeff2*B5;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 249.896909254999;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7_6;
    % Computing T2k9 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k9 = coeff1*T2k5 + coeff2*B7;
    output = T2k9;
end

