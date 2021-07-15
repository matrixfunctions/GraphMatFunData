n=MATSIZE;
A0=zeros(n,n);
for i=1:n
    for j=1:n
        if (i ~= j)
            A0(i,j)=1/sqrt(abs((i-j)));
        end
    end
end
A0(n,1) = A0(n,1) + 0.0001; % Break symmetry to avoid special case code
A0_basic=A0/norm(A0);


fprintf("System: %s  Computer: %s\n",version,computer);
addpath('/tmp');
fprintf("Matrix size: %d x %d \n", n , n);
version -blas ; blas=ans;
fprintf("BLAS version: %s\n\n",blas);

nof_samples=10; tv=zeros(nof_samples,1);

% Convenience functions
expm_matlab=@(x) expm(x);
expmpoly_matlab=@(x) expmpol(x);

if (~(exist("expmpol")>0))
    fprintf("Function expmpol() not in PATH. Trying to download from\n");
    fprintf("    http://personales.upv.es/~jorsasma/software/expmpol.m\n");
    fprintf("and save in current directory.... ");
    try
        websave("expmpol.m","http://personales.upv.es/~jorsasma/software/expmpol.m");
        fprintf("done.\n\n");
    catch
        fprintf("failed.\n")
        fprintf("Unable to download. expmpoly simulation set to identity.\n\n");
        expmpoly_matlab=@(x) x; % Dummy
    end
end



%%% START REPEATED CODE
% Code to run: NAME *******
fprintf("%20.20s: ","NAME");
A=MATNORM*A0_basic;

for k=1:nof_samples
    tic;
    NAME(A);
    tv(k)=toc;
    pause(0.5);
end
fprintf("%.6f ± %.3f\n",median(tv),std(tv));
fprintf("                       timings:");
fprintf("%.3f ",tv);
fprintf("\n");
pause(3)

%%% END REPEATED CODE


fprintf("Done!\n ");
exit
