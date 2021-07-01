using GraphMatFun, LinearAlgebra,GenericSVD

function graph_generalized_denman_beavers()

    # Construct a Degopt and change :mult to :ldiv
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

function normalize_denman_beavers(graph)
    (HA,HB,y)=get_coeffs_denman_beavers(graph);

    for k=1:size(HA,1)
        alpha=maximum(abs.(HB[k,:]));
        HA[k,:] /= alpha
        HB[k,:] /= alpha
        beta=maximum(abs.(HA[k,:]));
        HA[k,:] /= beta;
        if (k+2 <= size(HA,2))
            HA[:,k+2] /= beta;
            HB[:,k+2] /= beta;
        end
        y[k+2] /= beta;
    end

    degopt=Degopt(HA,HB,y);
    (graph,cref)=graph_degopt(degopt);
    # Replace :mult with ldiv
    for k in keys(graph.operations)
        if (graph.operations[k] == :mult)
            graph.operations[k] = :ldiv
        end
    end

    return (graph,cref)


end


include("newsimulationtools.jl");

statelist=[];

f=s-> sqrt(s+1);


## Standard Denman-Beavers with shift *********************
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





## Generalized Denman-Beavers (with shift) *********************
#T=ComplexF64;
state=State(f,complex(T));
state.params[:rho]=0.5;
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
                  linlsqr=:real_svd,droptol=1e-8,
                  stoptol=1e-16,errtype=:relerr,maxit=5)
showerr(state,showmeta=true)
opt_gauss_newton!(graph,f,discr,logger=1,cref=cref,
                  linlsqr=:real_svd,droptol=1e-9,
                  stoptol=1e-16,errtype=:relerr,maxit=5)
showerr(state,showmeta=true)
opt_gauss_newton!(graph,f,discr,logger=1,cref=cref,
                  linlsqr=:real_svd,droptol=1e-11,
                  stoptol=1e-16,errtype=:relerr,maxit=5)
showerr(state,showmeta=true)
opt_gauss_newton!(graph,f,discr,logger=1,cref=cref,
                  linlsqr=:real_svd,droptol=1e-12,
                  stoptol=1e-16,errtype=:relerr,maxit=5)
showerr(state,showmeta=true)
opt_gauss_newton!(graph,f,discr,logger=1,cref=cref,
                  linlsqr=:real_svd,droptol=1e-13,
                  stoptol=1e-16,errtype=:relerr,maxit=5)
showerr(state,showmeta=true)
opt_gauss_newton!(graph,f,discr,logger=1,cref=cref,
                  linlsqr=:real_svd,droptol=1e-14,
                  stoptol=1e-16,errtype=:relerr,maxit=10)
count(values(graph.operations) .== :ldiv)
push!(statelist,deepcopy(state))

for state in deepcopy(statelist)
    #state.graph=Compgraph(Float64,state.graph)
    showerr(state,showmeta=true)
end




## PS  *********************

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

## PS with degopt freedom optimization   *********************


x=get_coeffs(graph1,cref1);  x=x .+ 0.00001*rand(size(x)); set_coeffs!(graph1,x,cref1);

state1.params[:graphname]="ps_opt";
discr=get_disc_discr(state1);
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                  linlsqr=:real_svd,droptol=1e-5,γ0=0.5,
                  stoptol=1e-14,errtype=:relerr,maxit=10)
showerr(state1,showmeta=true)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                  linlsqr=:real_svd,droptol=1e-6,γ0=0.3,
                  stoptol=1e-14,errtype=:relerr,maxit=10)
showerr(state1,showmeta=true)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                         linlsqr=:real_svd,droptol=1e-7,γ0=0.1,
                         stoptol=1e-14,errtype=:relerr,maxit=10)
showerr(state1,showmeta=true)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                linlsqr=:real_svd,droptol=1e-7,
                                stoptol=1e-14,errtype=:relerr,maxit=20)
showerr(state1,showmeta=true)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                linlsqr=:real_svd,droptol=1e-8,γ0=0.1,
                                stoptol=1e-14,errtype=:relerr,maxit=20)
showerr(state1,showmeta=true)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                linlsqr=:real_svd,droptol=1e-8,γ0=0.5,
                                stoptol=1e-14,errtype=:relerr,maxit=20)
showerr(state1,showmeta=true)
#opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
#                  linlsqr=:real_svd,droptol=1e-9,γ0=0.5,
#                  stoptol=1e-14,errtype=:relerr,maxit=10)
#showerr(state1,showmeta=true)
push!(statelist,deepcopy(state1))

using Random
Random.seed!(0);
A=randn(100,100);
A=0.5*A/opnorm(A);

for state in deepcopy(statelist)
    state.graph=Compgraph(Float64,state.graph)
    showerr(state,showmeta=true)

    graph=state.graph;
    if (:ldiv in values(graph.operations))

        try
            (graph,_)=normalize_denman_beavers(graph)
            (HA,HB,y)=get_coeffs_denman_beavers(graph)

        catch
            HA=[];
            HB=[];
            y=[];
        end

    else
        (HA,HB,y)=get_degopt_coeffs(graph)
    end

    #@show latexify([HA HB])
    #@show latexify(y')

end


@show norm(eval_graph(graph,A)-sqrt(A+I))
@show norm(eval_graph(graph1,A)-sqrt(A+I))


saveall=false;
if (saveall)
    for state in deepcopy(statelist)
        graph=state.graph;
        if (:ldiv in values(graph.operations))

            try
                (graph,_)=normalize_denman_beavers(graph)
                (HA,HB,y)=get_coeffs_denman_beavers(graph)

            catch
            end
        end
        graphname=state.params[:graphname]
        basename=joinpath("..","..","data","sqrt_plus_one",graphname);
        @show "$(basename).cgr"
        export_compgraph(graph,"$(basename).cgr")
    end
end
