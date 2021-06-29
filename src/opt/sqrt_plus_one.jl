
using GraphMatFun, LinearAlgebra,GenericSVD
include("newsimulationtools.jl");

f=s-> sqrt(s+1);
(graph,cref)=graph_denman_beavers(2);
rename_node!(graph,:A,:A_shift,cref);
add_lincomb!(graph,:A_shift,1.0,:A,1.0,:I);
push!(cref,(:A_shift,1))
push!(cref,(:A_shift,2))

x=0.5;
f(x)-eval_graph(graph,x);
graph.outputs
o=graph.outputs[end];
empty!(graph.outputs)
add_output!(graph,o);
compress_graph!(graph,cref)



#T=ComplexF64;
T=BigFloat;
graph=Compgraph(T,graph);
state=State(f,complex(T));
state.params[:rho]=0.5;
state.params[:target_n]=2000;
state.params[:n]=12;
state.params[:graphname]="denman_beavers";
state.graph=graph;
showerr(state,showmeta=true)

discr=get_disc_discr(state);
opt_gauss_newton!(graph,f,discr,logger=1,cref=cref,
                  linlsqr=:real_svd,droptol=1e-12,
                  stoptol=1e-14,errtype=:relerr)
showerr(state,showmeta=true)

count(values(graph.operations) .== :ldiv)



(graph1,cref1)=graph_monomial(1:10.0)
#(graph1,cref1)=graph_monomial(1:6.0)
graph1=Compgraph(ComplexF64,graph1);
opt_linear_fit!(graph1,f,discr,cref1)
coeffs=get_coeffs(graph1,cref1);
(graph1,cref1)=graph_ps_degopt(real.(coeffs));

count(values(graph1.operations) .== :mult)
state1=deepcopy(state);
state1.graph=graph1;
state1.params[:n]=200;
state1.params[:graphname]="ps";
showerr(state1,showmeta=true)


state1.params[:graphname]="ps+opt";
discr=get_disc_discr(state1);
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                  linlsqr=:real_svd,droptol=1e-5,γ0=0.5,
                  stoptol=1e-14,errtype=:relerr,maxit=10)
graph1x=deepcopy(graph1);
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                  linlsqr=:real_svd,droptol=1e-6,γ0=0.3,
                  stoptol=1e-14,errtype=:relerr,maxit=10)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                         linlsqr=:real_svd,droptol=1e-7,γ0=0.1,
                         stoptol=1e-14,errtype=:relerr,maxit=100)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                linlsqr=:real_svd,droptol=1e-7,γ0=0.5,
                                stoptol=1e-14,errtype=:relerr,maxit=20)

opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                         linlsqr=:real_svd,droptol=1e-8,γ0=0.1,
                         stoptol=1e-14,errtype=:relerr,maxit=100)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                linlsqr=:real_svd,droptol=1e-7,γ0=0.5,
                                stoptol=1e-14,errtype=:relerr,maxit=20)


showerr(state1,showmeta=true)
