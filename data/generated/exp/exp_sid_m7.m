function output = exp_sid_m7(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: T2k2 B2 T2k3 Ba8_3 Bb6_3 Ba7_3 Bb8_3 B3 Bb6_4 Bb7_4 T2k4 Ba8_4 Bb8_4 Ba7_4 B4 Ba8_5 Bb8_5 T2k5 Bb6_5 Bb7_5 Ba7_5 B5 Ba8_6 Bb7_6 Ba7_6 Bb6 T2k6 Bb8_6 B6 Ba8_7 Ba7 Bb7 B7 Ba8 Bb8_7 B8 T2k10
    % Computing T2k2 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.5029302610017967;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Ba8_3 with operation: lincomb
    coeff1 = 11.10689398085882;
    coeff2 = 2.991654767354374;
    Ba8_3 = coeff1*A + coeff2*B2;
    % Computing Bb6_3 with operation: lincomb
    coeff1 = 1.313681421698863e-6;
    coeff2 = 6.204734935438909e-8;
    Bb6_3 = coeff1*A + coeff2*B2;
    % Computing Ba7_3 with operation: lincomb
    coeff1 = 5.951585263506065;
    coeff2 = 0.4155284057336423;
    Ba7_3 = coeff1*A + coeff2*B2;
    % Computing Bb8_3 with operation: lincomb
    coeff1 = -0.000263823622233776;
    coeff2 = 0.008139086096860678;
    Bb8_3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = A * B2;
    % Computing Bb6_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 2.957106114715868e-9;
    Bb6_4 = coeff1*Bb6_3 + coeff2*B3;
    % Computing Bb7_4 with operation: lincomb
    coeff1 = 0.03306559506631931;
    coeff2 = 0.002630043177655382;
    Bb7_4 = coeff1*B2 + coeff2*B3;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.07705596948494946;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Ba8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.2857950268422422;
    Ba8_4 = coeff1*Ba8_3 + coeff2*B3;
    % Computing Bb8_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.001121744731945438;
    Bb8_4 = coeff1*Bb8_3 + coeff2*B3;
    % Computing Ba7_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.02479095151834799;
    Ba7_4 = coeff1*Ba7_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = A * B3;
    % Computing Ba8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.03005135891320298;
    Ba8_5 = coeff1*Ba8_4 + coeff2*B4;
    % Computing Bb8_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 9.027588625491207e-5;
    Bb8_5 = coeff1*Bb8_4 + coeff2*B4;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.004985549176118462;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb6_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.556371639324141e-10;
    Bb6_5 = coeff1*Bb6_4 + coeff2*B4;
    % Computing Bb7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0002100333647757715;
    Bb7_5 = coeff1*Bb7_4 + coeff2*B4;
    % Computing Ba7_5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.001283057135586989;
    Ba7_5 = coeff1*Ba7_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = A * B4;
    % Computing Ba8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.002742336655922557;
    Ba8_6 = coeff1*Ba8_5 + coeff2*B5;
    % Computing Bb7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.7537107416419e-5;
    Bb7_6 = coeff1*Bb7_5 + coeff2*B5;
    % Computing Ba7_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 3.501669195497238e-5;
    Ba7_6 = coeff1*Ba7_5 + coeff2*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.556371639324141e-11;
    Bb6 = coeff1*Bb6_5 + coeff2*B5;
    % Computing T2k6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 6.263526066651383e-5;
    T2k6 = coeff1*T2k5 + coeff2*B5;
    % Computing Bb8_6 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 8.572383602707347e-6;
    Bb8_6 = coeff1*Bb8_5 + coeff2*B5;
    % Computing B6 with operation: mult
    B6 = B5 * Bb6;
    % Computing Ba8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 61.75954247606858;
    Ba8_7 = coeff1*Ba8_6 + coeff2*B6;
    % Computing Ba7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Ba7 = coeff1*Ba7_6 + coeff2*B6;
    % Computing Bb7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Bb7 = coeff1*Bb7_6 + coeff2*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Ba8 = coeff1*Ba8_7 + coeff2*B7;
    % Computing Bb8_7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    Bb8_7 = coeff1*Bb8_6 + coeff2*B6;
    % Computing B8 with operation: mult
    B8 = Ba8 * Bb8_7;
    % Computing T2k10 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    T2k10 = coeff1*T2k6 + coeff2*B8;
    output = T2k10;
end

