function output = exp_sastre_m4_opt_rho0_69(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Bb4_3 Ba3_2 Ba3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 T2k5 Bb5 B5 T2k7
    % Computing Bb2 with operation: lincomb
    coeff1 = -0.26484760452432327;
    coeff2 = 0.9903047529247974;
    Bb2 = coeff1*I + coeff2*A;
    % Computing Ba4_2 with operation: lincomb
    coeff1 = 0.0002727973329217783;
    coeff2 = 0.05658254301026135;
    Ba4_2 = coeff1*I + coeff2*A;
    % Computing Bb4_2 with operation: lincomb
    coeff1 = 0.0045397184361004065;
    coeff2 = 0.0026543192111867;
    Bb4_2 = coeff1*I + coeff2*A;
    % Computing Ba2 with operation: lincomb
    coeff1 = -0.26484760452432327;
    coeff2 = 0.9903047529247974;
    Ba2 = coeff1*I + coeff2*A;
    % Computing B2 with operation: mult
    B2 = Ba2 * Bb2;
    % Computing Ba4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = -0.31616784713314605;
    Ba4_3 = coeff1*Ba4_2 + coeff2*B2;
    % Computing Bb4_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.0004235929351940947;
    Bb4_3 = coeff1*Bb4_2 + coeff2*B2;
    % Computing Ba3_2 with operation: lincomb
    coeff1 = -0.3183560628971706;
    coeff2 = 0.9451759356891595;
    Ba3_2 = coeff1*I + coeff2*A;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.08341693362546203;
    Ba3 = coeff1*Ba3_2 + coeff2*B2;
    % Computing Ba5_2 with operation: lincomb
    coeff1 = 0.10536695630497932;
    coeff2 = 1.3909408033494417;
    Ba5_2 = coeff1*I + coeff2*A;
    % Computing Ba5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.3494801360488722;
    Ba5_3 = coeff1*Ba5_2 + coeff2*B2;
    % Computing Bb5_2 with operation: lincomb
    coeff1 = 0.14340598505310592;
    coeff2 = 0.18419588627720188;
    Bb5_2 = coeff1*I + coeff2*A;
    % Computing Bb5_3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.06331721212120305;
    Bb5_3 = coeff1*Bb5_2 + coeff2*B2;
    % Computing T2k2 with operation: lincomb
    coeff1 = 0.9584219043602061;
    coeff2 = 1.0556647430407748;
    T2k2 = coeff1*I + coeff2*A;
    % Computing T2k3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.4759594481148073;
    T2k3 = coeff1*T2k2 + coeff2*B2;
    % Computing Bb3_2 with operation: lincomb
    coeff1 = 0.041562405017987775;
    coeff2 = -0.30378066944691157;
    Bb3_2 = coeff1*I + coeff2*A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.949941947908575;
    Bb3 = coeff1*Bb3_2 + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.9736101476608433;
    Ba4 = coeff1*Ba4_3 + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 4.468866390298052e-5;
    Bb4 = coeff1*Bb4_3 + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Ba5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.029870845474797018;
    Ba5_4 = coeff1*Ba5_3 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.4576696798781078;
    Ba5 = coeff1*Ba5_4 + coeff2*B4;
    % Computing T2k4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.11560289591544066;
    T2k4 = coeff1*T2k3 + coeff2*B3;
    % Computing Bb5_4 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.016928981566343947;
    Bb5_4 = coeff1*Bb5_3 + coeff2*B3;
    % Computing T2k5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 5.152589478327363;
    T2k5 = coeff1*T2k4 + coeff2*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.6928011352478407;
    Bb5 = coeff1*Bb5_4 + coeff2*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing T2k7 with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 0.7070319120970353;
    T2k7 = coeff1*T2k5 + coeff2*B5;
    output = T2k7;
end

