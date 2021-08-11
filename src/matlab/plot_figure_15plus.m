% Generates the figure 5. in the manuscript https://arxiv.org/abs/2107.12198

% Download from https://github.com/matrixfunctions/GraphMatFunData/tree/main/data/exp
degopt1=read_degopt('exp_sastre_m4.cgr',true);
degopt2=read_degopt('exp_sastre_m4_opt_rho0_69.cgr',true);

xv=0:0.05:1.0;
xv(1) = xv(2)/10;
err1=vpa(zeros(size(xv)));
err2=vpa(zeros(size(xv)));

for i=1:length(xv)
    x=vpa(xv(i));
    err1(i)=(eval_degopt(degopt1,x)-exp(x))/exp(x);
    err2(i)=(eval_degopt(degopt2,x)-exp(x))/exp(x);
end

semilogy(xv,double(abs(err1)),xv,double(abs(err2)));

legend('15+', 'opt');