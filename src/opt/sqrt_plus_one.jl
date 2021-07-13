using GraphMatFun, LinearAlgebra,GenericSVD, Random
include("simulationtools.jl");

# Julia's base can change the random seed generators.
# We hardcode it to make the results reproducible over
# Julia versions
function hardcoded_rand(t)

    x=[  0.8872174233685589;  0.8377413644286786;  0.5608614162945871;  0.7142032679122071;  0.8217070528857874;  0.40399577722972235;  0.8659919887404998;  0.08029198288895634;  0.10052875305338604;  0.40627896679851794;  0.8342948272930204;  0.3820649650758666;  0.016757502407812153;  0.20137662172350335;  0.2555463313583042;  0.19248011642467033;  0.011618055692006912;  0.42854938782243956;  0.21254339113365306;  0.8654170264281245;  0.26087081576693616;  0.0330965543051257;  0.15171093545767034;  0.541356766187003;  0.7025653111402541;  0.3832283012090676;  0.24532744474494972;  0.8773921298355143;  0.28729058543087427;  0.02645861903949709;  0.10877066841853389;  0.25559799376084;  0.9817655834384857;  0.09350052884595439;  0.5129871676353079;  0.9371434490087653;  0.8916549976510515;  0.9426774471099136;  0.45772957349208;  0.9332990364480496;  0.8382243800948085;  0.40812659546896835;  0.653899520291725;  0.6110199290298604;  0.19036516891430988;  0.26967359818868164;  0.7360989887783765; 0.296931260360374];
    x=x[1:t[1]];
    return x;

end

# Construct Denman-Beavers iteration for the modified
# degree-optimal form.
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
# Get all the coefficient valeus in the mod degopt stored in graph
function get_coeffs_mod_degopt(graph)
    graph=deepcopy(graph);
    for k in keys(graph.operations)
        if (graph.operations[k] == :ldiv)
            graph.operations[k] = :mult
        end
    end
    return get_degopt_coeffs(graph);
end
# Transform all rational functions such with a normalization
# of the maximum coefficient value in order to avoid
# large/small numbers
function normalize_mod_degopt!(graph)
    (HA,HB,y)=get_coeffs_mod_degopt(graph);

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
println("Optimizing Denman-Beavers (modified degree-optimal form) *********************");
T=BigFloat;
state=State(f,complex(T));
state.params[:rho]=0.5;
state.params[:target_n]=2000;
state.params[:n]=1000;
state.params[:graphname]="denman_beavers";
state.graph=graph;
showerr(state,showmeta=true)
push!(statelist,deepcopy(state))

dd=count(values(graph.operations) .== :ldiv)

println("Using $dd inverse multiplications");

x=get_coeffs(graph,cref);  x=x .+ 0.00001*hardcoded_rand(size(x)); set_coeffs!(graph,x,cref);

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




println("Optimizing PS  *********************");

# Computing the least squares match of coefficients on the disk
#(graph1,cref1)=graph_monomial(1:10.0) # If we want 4 multiplications
(graph1,cref1)=graph_monomial(1:13.0)
graph1=Compgraph(ComplexF64,graph1);
opt_linear_fit!(graph1,f,discr,cref1)
coeffs=get_coeffs(graph1,cref1);
(graph1,cref1)=graph_ps_degopt(real.(coeffs));

mm=count(values(graph1.operations) .== :mult)
println("Using $mm multiplications");
state1=deepcopy(state);
state1.graph=graph1;
state1.params[:n]=200;
state1.params[:graphname]="ps_m$(mm)";
showerr(state1,showmeta=true)
push!(statelist,deepcopy(state1))

## PS with degopt freedom optimization   *********************


x=get_coeffs(graph1,cref1);  x=x .+ 0.00001*hardcoded_rand(size(x,1)); set_coeffs!(graph1,x,cref1);

state1.params[:graphname]="ps_m$(mm)_opt";
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
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                       linlsqr=:real_svd,droptol=1e-8,γ0=1.0,
                                       stoptol=1e-14,errtype=:relerr,maxit=5)

showerr(state1,showmeta=true)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                linlsqr=:real_svd,droptol=1e-9,γ0=0.5,
                  stoptol=1e-14,errtype=:relerr,maxit=20)
showerr(state1,showmeta=true)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                  linlsqr=:real_svd,droptol=1e-9,γ0=1.0,
                  stoptol=1e-14,errtype=:relerr,maxit=5)
showerr(state1,showmeta=true)

opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                                       linlsqr=:real_svd,droptol=1e-10,γ0=0.5,
                         stoptol=1e-16,errtype=:relerr,maxit=30)
showerr(state1,showmeta=true)
opt_gauss_newton!(graph1,f,discr,logger=1,cref=cref1,
                  linlsqr=:real_svd,droptol=1e-10,γ0=1.0,
                  stoptol=1e-14,errtype=:relerr,maxit=5)



push!(statelist,deepcopy(state1))

Random.seed!(0);
A=randn(100,100);
A=0.5*A/opnorm(A);

for state in deepcopy(statelist)
    state.graph=Compgraph(Float64,state.graph)
    showerr(state,showmeta=true)

    graph=state.graph;
    if (:ldiv in values(graph.operations))

        try
            (graph,_)=normalize_mod_degopt!(graph)
            (HA,HB,y)=get_coeffs_mod_degopt(graph)

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
                (graph,_)=normalize_mod_degopt!(graph)
                (HA,HB,y)=get_coeffs_mod_degopt(graph)

            catch
            end
        end
        graphname=state.params[:graphname]
        basename=joinpath("..","..","data","sqrt_plus_one",graphname);
        @show "$(basename).cgr"
        export_compgraph(graph,"$(basename).cgr")
    end
end
