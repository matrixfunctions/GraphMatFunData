using GraphMatFun, LinearAlgebra,GenericSVD

function graph_generalized_denman_beavers()
    graph=Compgraph(Float64);

    HA=[1.0  1.0   0.0  0.0  0.0
        0.5  0.0   0.5  0.0  0.0
        1.0  0.5   0.0  0.0  0.0
        0.25  0.0  0.25 0.0  0.5];

    HB=[1.0 0.0 0.0 0.0 0.0
        1.0 0.0 0.0 0.0 0.0
        1.0 0.0 0.0 0.0 0.0
        1.0 0.0 0.0 0.0 0.0];
    y=[1/4; 1/8;0;1/4;0;1/2];


    degopt=Degopt(HA,HB,y);
    (graph,cref)=graph_degopt(degopt);
    # Replace :mult with ldiv
    for k in keys(graph.operations)
        if (graph.operations[k] == :mult)
            graph.operations[k] = :ldiv
        end
    end

#
#    cref=[];
#    cref_tmp=add_sum!(graph,:B1a,[1; 1.0], [:I;:A]);;
#    push!(cref,cref_tmp...);
#    cref_tmp=add_sum!(graph,:B1b,[1; 0.0], [:I;:A]);;
#    push!(cref,cref_tmp...);
#    add_ldiv!(graph,:B1,:B1a,:B1b);
#
#
#    cref_tmp=add_sum!(graph,:B2a,[0.5; 0; 0.5], [:I;:A;:B1]);;
#    push!(cref,cref_tmp...);
#    cref_tmp=add_sum!(graph,:B2b,[1; 0;0.0], [:I;:A;:B1]);;
#    push!(cref,cref_tmp...);
#    add_ldiv!(graph,:B2,:B2a,:B2b);
#
#
#    cref_tmp=add_sum!(graph,:B3a,[1; 0.5;0;0], [:I;:A;:B1;:B2]);;
#    push!(cref,cref_tmp...);
#    cref_tmp=add_sum!(graph,:B3b,[1; 0;0;0.0], [:I;:A;:B1;:B2]);;
#    push!(cref,cref_tmp...);
#    add_ldiv!(graph,:B3,:B3a,:B3b);
#
#
#    cref_tmp=add_sum!(graph,:B4a,[1/4; 0; 1/4; 0; 1/2], [:I;:A;:B1;:B2;:B3]);;
#    push!(cref,cref_tmp...);
#    cref_tmp=add_sum!(graph,:B4b,[1; 0;0;0;0.0], [:I;:A;:B1;:B2;:B3]);;
#    push!(cref,cref_tmp...);
#    add_ldiv!(graph,:B4,:B4a,:B4b);
#
#    cref_tmp=add_sum!(graph,:y,[1/4; 1/8;0;1/4;0;1/2],
#                      [:I;:A;:B1;:B2;:B3;:B4])
#    push!(cref,cref_tmp...);
#    add_output!(graph,:y);
#
    return (graph,cref)
end
function get_coeffs_denman_beavers(graph)
    graph=deepcopy(graph);
    for k in keys(graph.operations)
        if (graph.operations[k] == :ldiv)
            graph.operations[k] = :mult
        end
    end

    return get_degopt_coeffs(graph);

end


include("newsimulationtools.jl");

statelist=[];

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
state.params[:rho]=0.4;
state.params[:target_n]=2000;
state.params[:n]=1000;
state.params[:graphname]="denman_beavers";
state.graph=graph;
showerr(state,showmeta=true)
push!(statelist,deepcopy(state))


x=get_coeffs(graph,cref);  x=x .+ 0.00001*rand(size(x)); set_coeffs!(graph,x,cref);

(graph,cref)= graph_generalized_denman_beavers();
state.graph=graph;
state.params[:graphname]="denman_beavers_opt";
discr=get_disc_discr(state);

opt_gauss_newton!(graph,f,discr,logger=1,cref=cref,
                  linlsqr=:real_svd,droptol=1e-10,
                  stoptol=1e-14,errtype=:relerr,maxit=10)
showerr(state,showmeta=true)
opt_gauss_newton!(graph,f,discr,logger=1,cref=cref,
                  linlsqr=:real_svd,droptol=1e-14,
                  stoptol=1e-14,errtype=:relerr,maxit=10)
count(values(graph.operations) .== :ldiv)
push!(statelist,deepcopy(state))

for state in deepcopy(statelist)
    #state.graph=Compgraph(Float64,state.graph)
    showerr(state,showmeta=true)
end

get_coeffs_denman_beavers(graph,cref)
asd
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
push!(statelist,deepcopy(state1))


state1.params[:graphname]="ps+opt0";
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
push!(statelist,deepcopy(state1))
state1.params[:graphname]="ps+opt1";
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                         linlsqr=:real_svd,droptol=1e-8,γ0=0.1,
                         stoptol=1e-14,errtype=:relerr,maxit=100)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                linlsqr=:real_svd,droptol=1e-8,γ0=0.5,
                                stoptol=1e-14,errtype=:relerr,maxit=20)

opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                       linlsqr=:real_svd,droptol=1e-9,γ0=0.1,
                                       stoptol=1e-14,errtype=:relerr,maxit=20)
push!(statelist,deepcopy(state1))
state1.params[:graphname]="ps+opt2";

opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                       linlsqr=:real_svd,droptol=1e-9,γ0=0.5,
                                       stoptol=1e-14,errtype=:relerr,maxit=20)

opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                       linlsqr=:real_svd,droptol=1e-10,γ0=0.1,
                                       stoptol=1e-14,errtype=:relerr,maxit=20)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                       linlsqr=:real_svd,droptol=1e-10,γ0=0.5,
                                       stoptol=1e-14,errtype=:relerr,maxit=30)

push!(statelist,deepcopy(state1))
state1.params[:graphname]="ps+opt3";

opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                  linlsqr=:real_svd,droptol=3e-11,γ0=0.01,
                  stoptol=1e-14,errtype=:relerr,maxit=40)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                  linlsqr=:real_svd,droptol=3e-11,γ0=0.5,
                  stoptol=1e-14,errtype=:relerr,maxit=40)


opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                       linlsqr=:real_svd,droptol=5e-11,γ0=0.1,
                                       stoptol=1e-14,errtype=:relerr,maxit=20)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                       linlsqr=:real_svd,droptol=1e-11,γ0=0.5,
                                       stoptol=1e-14,errtype=:relerr,maxit=30)
push!(statelist,deepcopy(state1))
state1.params[:graphname]="ps+opt4";

#opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
#                                linlsqr=:real_svd,droptol=1e-9,γ0=0.5,
#                                stoptol=1e-14,errtype=:relerr,maxit=60)
#
push!(statelist,deepcopy(state1))
showerr(state1,showmeta=true)

using Random
Random.seed!(0);
A=randn(100,100);
A=0.3*A/opnorm(A);

for state in deepcopy(statelist)
    state.graph=Compgraph(Float64,state.graph)
    showerr(state,showmeta=true)
end


@show norm(eval_graph(graph,A)-sqrt(A+I))
@show norm(eval_graph(graph1,A)-sqrt(A+I))
