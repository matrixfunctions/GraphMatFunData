function output=denman_beavers(A)
    n=size(A,1);
    I=eye(n,n);
    % Computation order: A_shift Xinv0 Y1 X1 Yinv1 X2 Xinv1 Y2 Yinv2 X3
    % Computing A_shift with operation: lincomb
    coeff1=1.0;
    coeff2=1.0;
    A_shift= coeff1*A+coeff2*I;
    % Computing Xinv0 with operation: ldiv
    Xinv0=A_shift \ I;
    % Computing Y1 with operation: lincomb
    coeff1=0.5;
    coeff2=0.5;
    Y1= coeff1*I+coeff2*Xinv0;
    % Computing X1 with operation: lincomb
    coeff1=0.5;
    coeff2=0.5;
    X1= coeff1*A_shift+coeff2*I;
    % Computing Yinv1 with operation: ldiv
    Yinv1=Y1 \ I;
    % Computing X2 with operation: lincomb
    coeff1=0.5;
    coeff2=0.5;
    X2= coeff1*X1+coeff2*Yinv1;
    % Computing Xinv1 with operation: ldiv
    Xinv1=X1 \ I;
    % Computing Y2 with operation: lincomb
    coeff1=0.5;
    coeff2=0.5;
    Y2= coeff1*Y1+coeff2*Xinv1;
    % Computing Yinv2 with operation: ldiv
    Yinv2=Y2 \ I;
    % Computing X3 with operation: lincomb
    coeff1=0.5;
    coeff2=0.5;
    X3= coeff1*X2+coeff2*Yinv2;
    output=X3;
end

