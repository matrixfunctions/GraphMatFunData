using GenericLinearAlgebra, DelimitedFiles,Random
# State of the
mutable struct State
    f   # function
    discr
    eltype
    #
    graph
    cref
    #
    params # Dict with keywords
end
function State(f,eltype)
    return State(f,[],eltype,Nothing,Nothing,Dict{Symbol,Any}(:opt_kwargs => Dict()));
end

function scale_and_square!(state)
    degopt=Degopt(state.graph)
    scale!(degopt,1/2)
    square!(degopt)
    (g,c)=graph_degopt(degopt)
    state.graph=g
    state.cref=c
end


function get_disc_discr(state;n=state.params[:n])
    rho=state.params[:rho]
    nn=round(Int,n);
    discr=rho*exp.(1im*range(0,2*pi,length=nn)[1:end-1]);
    discr=convert.(state.eltype,discr);
    return discr;
end
function get_taylor(eltype,n)
    c=1 ./factorial.(convert.(big(real(eltype)),0:(n-1)));
    c=convert.(eltype,c);
end


# Corresponds to a call to opt_gaussnewton!
struct OptSimulation
    kwargs
    cref
    special_case
end

function OptSimulation(s::State;special_case=0)
    return OptSimulation(s.params[:opt_kwargs],s.cref,special_case);
end

# Levenberg-Marquard simulations using LsqFit.curve_fit
struct LMSimulation
    kwargs
end

struct NoOp
end

function LMSimulation(s::State)
    kwargs=Dict();
    if (haskey(s.params,:curve_fit_kwargs))
        kwargs=s.params[:curve_fit_kwargs]
    end
    return LMSimulation(kwargs)
end


# Corresponds to changing a parameter
struct ModifyParam
    param
    factor
end

# Corresponds to a change of all variables
struct KickIt
    factor::Number
    mode
end
function KickIt(state::State,factor)

    if haskey(state.params,:kickit_mode)
        mode=state.params[:kickit_mode]
    else
        mode=0;
    end
    @show mode
    return KickIt(factor,mode);

end
struct ZeroNormalize
end

struct ShowInfo
    theta_guess;
end

struct LineSearch
end

struct SetSmallToZero
    k
end
struct PurgeZeroOpt
    k
end

struct DegoptNormalize
end

# Initialize State with a graph based on number of multiplications.
function init_state_mult(f,rho,n,mult; eltype=Complex{BigFloat})
    init_state_mult(nothing,f,rho,n,mult; eltype=eltype);
end

function init_state_mult(graphname,f,rho,n,mult; eltype=Complex{BigFloat})
    if (f != exp)
        error("Only exp implemented");
    end
    state=State(f,eltype);
    state.params[:n]=n;
    state.params[:target_n]=n;
    state.params[:rho]=rho;
    state.params[:target_n]=n;

   state.discr=get_disc_discr(state);


    init_state_mult!(state,graphname,mult);
    return state;
end

function init_state_mult!(state,graphname,mult;showmeta=false)

    if (!isnothing(graphname))
        print("Generating $graphname ");
        if (graphname == :sid)
            (graph,cref)=graph_sid_exp(mult)
        elseif (graphname == :bbc)
            (graph,cref)=graph_bbc_exp(m)
        elseif (graphname == :sastre)
            method=:y1s;
            if (m==6)
                method=:h2m
            elseif (m==8)
                method=:z1ps
            end
            (graph,cref)=graph_sastre_exp(m,method)
        end

        if (graphname in [:ps,:mono])

            taylorcoeffs=get_taylor(state.eltype,200);

            # Automatically detect degree
            graph_mult = 0;
            deg0 = m-1;
            graph=nothing; oldgraph=nothing;
            cref=nothing; oldcref=nothing;
            while graph_mult<mult+1
                deg0 += 1;
                c=taylorcoeffs[1:(deg0+1)];
                oldgraph=graph;
                oldcref=cref;
                if (graphname == :ps)
                    (graph,cref)=graph_ps_degopt(c)
                elseif (graphname == :mono)
                    (graph,cref)=graph_monomial_degopt(c)
                end
                graph_mult=count(values(graph.operations) .==:mult);
            end
            # We go one step too far in the automatic determination so
            # use the old graph
            graph=oldgraph;
            cref=oldcref;

            println(" degree $(deg0-1)");
        end

        state.params[:graphname]="$(graphname)_m$(mult)"
        graph = Compgraph(state.eltype,graph);
        state.graph=graph;
        state.cref=cref;
        if (showmeta)
            showerr(state,showmeta=true)
        end

    end
    return state;
end
function init_state_file!(state,graphname,filename;showmeta=true,scale_and_square=false);

    if (filename isa String)
        state.graph=import_compgraph(filename);
    else
        state.graph=filename;
    end

    mult=count(values(state.graph.operations) .==:mult);
    if (scale_and_square)
        mult += 1;
    end
    state.params[:graphname]="$(graphname)_m$(mult)";
    if (scale_and_square)
        scale_and_square!(state)
    else
        (_,cref)=graph_degopt(Degopt(state.graph))
        state.cref=cref;
    end


    if (showmeta)
        println("Generating $graphname based on $filename");
        showerr(state,showmeta=true)
    end

    return state;
end

# initcommand(control::String,state)  --> command object
function initcommand(control,state)
    if (control == "s")
        return OptSimulation(state);
    elseif (control == "S")
        return OptSimulation(state;special_case=1);
    elseif (control == "2")
        return OptSimulation(state;special_case=2);
    elseif (control == "3")
        return OptSimulation(state;special_case=3);
    elseif (control == "4")
        return OptSimulation(state;special_case=4);
    elseif (control == "5")
        return OptSimulation(state;special_case=5);
    elseif (control == "l")
        return LMSimulation(state);
    elseif (control == "L")
        return LineSearch();
    elseif (control == "g")
        return ModifyParam(:γ0,0.8);
    elseif (control == "G")
        return ModifyParam(:γ0,1/0.8);
    elseif (control=="n")
        return ModifyParam(:n,1/2);
    elseif (control=="Z")
        return SetSmallToZero(1);
    elseif (control=="N")
        return ModifyParam(:n,2);
    elseif (control=="D")
        return ModifyParam(:droptol,2);
    elseif (control=="d")
        return ModifyParam(:droptol,1/2);
    elseif (control=="k")
        return KickIt(state,1e-9);
    elseif (control=="K")
        return KickIt(state,1e-6);
    elseif (control=="z")
        return ZeroNormalize();
    elseif (control=="p") # Purge crefs
        return PurgeZeroOpt(1);
    elseif (control=="0") #
        return DegoptNormalize();
    elseif (control=="t")
        return ShowInfo(Nothing);
    elseif (control=="T")
        return ShowInfo(0.3);
    else
        return NoOp();
    end
end

# runcommand(commandobj, state) -> State
function runcommand(s::OptSimulation,state)
    # Run Gauss-Newton with params in s.params[:opt_kwargs]
    f=state.f;
    graph=deepcopy(state.graph);
    discr=state.discr;
    cref=deepcopy(s.cref);
    if (s.special_case == 1)
        discr=deepcopy(discr);
        push!(discr,zero(eltype(discr)));
    end
    if (s.special_case == 2)
        discr0=discr;
        discr=deepcopy(discr);
        append!(discr,discr0/2);
        append!(discr,0.1*discr0/maximum(abs.(discr0)));
    end
    if (s.special_case == 3)
        discr=discr/2;
    end
    if (s.special_case == 4)
        discr=discr/10;
    end
    if (s.special_case == 5)
        c=state.params[:cref_count]
        cref=[cref[mod(c,size(cref,1))+1]]
        state.params[:cref_count] += 1;
        println("count c:$c $(cref[1])");
    end


    opt_gauss_newton!(graph,f,discr;logger=1,
                      stoptol=1e-40,cref=cref,
                      s.kwargs...);

    ff=f.(discr);
    ftilde=eval_graph(graph,discr);
    errv=abs.((ff-ftilde)./ftilde);



    state=deepcopy(state);

    state.graph=graph # Set the output
    currentname=state.params[:graphname]
    if (!contains(currentname,"_opt"))
        state.params[:graphname]=currentname*"_opt";
    end

    return state;
end

function runcommand(s::ModifyParam,state)
    # Update a parameter some opt_kwargs treated separately
    state=deepcopy(state);
    param = s.param;
    factor = s.factor;
    if (param in [:droptol,:γ0])
        o=state.params[:opt_kwargs]
    else
        o=state.params;
    end
    T=typeof(o[param])
    o[param]=o[param]*factor;
    println("$param=$(o[param])");
    if (param == :n)
        state.discr=get_disc_discr(state);
    end

    return state;
end

function runcommand(s::LineSearch,state)

    graph0=deepcopy(state.graph);
    smallest_err=Inf;
    smallest_k=Inf;
    smallest_graph=deepcopy(state.graph);
    for (k,cref) in enumerate(state.cref)
        cref=state.cref[1];
        err0=showerr(state,output=false);


        v0=get_coeffs(state.graph,[cref])[1];

        vv=[v0];
        append!(vv,v0*(1 .+ 10 .^(-20.0:2:-1)))
        append!(vv,v0*(1 .- 10 .^(-20.0:2:-1)))
        #vv=v0*(1 .+ [0;])
        errv=zeros(typeof(err0),size(vv));
        for (i,v) in enumerate(vv)
            state.graph=deepcopy(graph0);
            set_coeffs!(state.graph,v,[cref]);
            errv[i]=showerr(state,output=false);
        end
        j=argmin(errv);
        state.graph=deepcopy(graph0);
        set_coeffs!(state.graph,vv[j],[cref]);
        if (errv[j] <= smallest_err)
            #@show errv
            smallest_err=errv[j];
            smallest_k=k;
            smallest_graph=deepcopy(state.graph);
        end
    end
    println("Modifying $(state.cref[smallest_k])");
    state.graph=smallest_graph;
    return state;

end
function runcommand(s::SetSmallToZero,state)
    num_elements=s.k
    state=deepcopy(state);
    graph=deepcopy(state.graph);
    #crefs=get_all_cref(graph);
    crefs=deepcopy(state.cref);
    vals=abs.(get_coeffs(graph,crefs));
    J=sortperm(vals)
    removed=Float64(norm(vals[J[1:num_elements]]));
    println("zeroing: $(removed)");
    set_coeffs!(graph,zero(vals[J[1:num_elements]]),
                crefs[J[1:num_elements]]);
    state.graph=graph;

    vals=abs.(get_coeffs(graph,crefs));
    J=findall(vals .== 0);

    return state;
end

function runcommand(s::PurgeZeroOpt,state)

    num_elements=s.k
    state=deepcopy(state);
    graph=deepcopy(state.graph);
    #crefs=get_all_cref(graph);
    crefs=deepcopy(state.cref);
    vals=abs.(get_coeffs(graph,crefs));
    J=findall(vals .== 0);
    if (size(J,1)>0)
        println("Purging cref $(crefs[J[1:s.k]])");
        setdiff!(state.cref,crefs[J[1:s.k]]);
    else
        println("Found nothing to purge");
    end
    return state;
end
function runcommand(s::DegoptNormalize,state)
    state=deepcopy(state);
    degopt=Degopt(state.graph);
    normalize!(degopt);
    (g,_)=graph_degopt(degopt);
    state.graph=g;
    return state;
end

function curve_fit_model(x,coeffs,graph,cref)
    p=size(x,1)/2;
    org_coeffs=get_coeffs(graph,cref);
    set_coeffs!(graph,coeffs,cref)
    z=eval_graph(graph,x);
    #set_coeffs!(graph,org_coeffs,cref);
    return [real.(z);imag.(z)];
end


function curve_fit_jacobian_model(x,coeffs,graph,cref)
    org_coeffs=get_coeffs(graph,cref);
    set_coeffs!(graph,coeffs,cref)
    J=eval_jac(graph,x,cref);
    #set_coeffs!(graph,org_coeffs,cref);
    return [real.(J);imag.(J)]
end

function runcommand(s::LMSimulation,state)

    discr=state.discr;
    T=real(state.eltype);
    xdata=complex(T).(discr);
    ydata_c=(state.f.(discr));
    ydata=T.([real.(ydata_c);imag(ydata_c)])

    graph=Compgraph(T,state.graph);
    cref=state.cref;
    p0=get_coeffs(graph,cref);


    model = (x,coeffs) -> curve_fit_model(x,coeffs,graph,cref)
    jacobian_model =
        (x,coeffs) -> curve_fit_jacobian_model(x,coeffs,graph,cref)

    fit = curve_fit(model, jacobian_model, xdata, ydata, p0,#show_trace=true,
                    x_tol=big(1e-20), g_tol=1e-20, maxIter=10,
                    s.kwargs...)

    state.graph=graph # Set the output
    @show typeof(graph)
    currentname=state.params[:graphname]
    if (!contains(currentname,"_opt"))
        state.params[:graphname]=currentname*"_opt";
    end
    return state;

end
function runcommand(s::NoOp,state)
    println("No op");
    return state;
end

function runcommand(s::KickIt,state)
    factor=s.factor;
    (x,y)=get_degopt_crefs(state.graph);
    random_state=FakeRandomState(0);
    #Random.seed!(0);
    if (s.mode == 0)
        cref=[x[end][1][1]];
        vals=get_coeffs(state.graph,cref);
        vals=vals.*(1 .+ fake_randn(random_state,size(vals,1)));
        set_coeffs!(state.graph,vals,cref);

    elseif (s.mode == 1)
        cref=[x[end][1][1]];
        vals=get_coeffs(state.graph,cref);
        vals=vals .+ fake_randn(random_state,size(vals,1));
        set_coeffs!(state.graph,vals,cref);
    elseif (s.mode == 2)

        cref=[x[1][2][1]];
        vals=get_coeffs(state.graph,cref);
        vals .+= 0.0001
        set_coeffs!(state.graph,vals,cref);
    elseif (s.mode == 3)
        cref=[x[end][2][1]];
        vals=get_coeffs(state.graph,cref);
        vals=vals .+ 0.0001*fake_randn(random_state,size(vals,1));
        set_coeffs!(state.graph,vals,cref);
    elseif s.mode ==4
        y_cref=get_degopt_crefs(state.graph)[2];
        # Kick the y-degopt vector
        Ha=[ -0.000599134   1.00256       0.0          0.0           0.0         0.0         0.0
 -0.000118639  -0.000778662   0.999709     0.0           0.0         0.0         0.0
 -1.09012f-5   -0.000246103  -0.00126751   0.999822      0.0         0.0         0.0
 -1.08452f-6   -3.92744f-5   -0.00017262  -0.00189144    0.999981    0.0         0.0
  5.84573f-7   -4.80433f-6   -2.1741f-5   -0.000177471  -0.00239572  0.999999    0.0
  1.25109f-7   -1.05495f-7   -2.48853f-6  -4.14857f-5    4.32646f-5  0.00116684  1.0
             ];
        Hb=[-0.000599134  1.00256   0.0        0.0         0.0          0.0         0.0
 -0.000838943  1.0011    0.0166752  0.0         0.0          0.0         0.0
 -0.0019134    1.00034   0.0377923  0.00100261  0.0          0.0         0.0
 -0.00242811   0.999187  0.0604852  0.00308486  6.53079f-5   0.0         0.0
  0.00116408   0.997916  0.0755767  0.0051393   0.000163854  6.79643f-8  0.0
 -8.09922f-8   0.997063  0.0954523  0.00867858  0.000495314  1.32541f-5  3.94245f-8]
        crand=[0.999999803981589
               1.0006179487544273
               0.4979281278966493
               0.16592558923333828
               0.03871151836039848
               0.0061833894649526416
               0.0006066317390421412
               3.0487323320460782e-5]
        (gtmp,cref)=graph_degopt(Degopt(Ha,Hb,crand));
        cvals=get_coeffs(gtmp,cref);
        set_coeffs!(state.graph,cvals,cref);
        #p=min(size(y_cref,1),size(crand,1));
        #set_coeffs!(state.graph,crand[1:p],y_cref);


    elseif (s.mode == 5)
        crefs=state.cref;
        vals=get_coeffs(state.graph,crefs);
        vals=vals .+ 0.0001*fake_randn(random_state,size(vals,1));
        set_coeffs!(state.graph,vals,crefs);
    else

        println("Unknown Kickit Mode");
    end

    return state;
end

function runcommand(s::ZeroNormalize,state)
    graph=state.graph;
    print("Diff before:"*string(eval_graph(graph,0)-state.f(0)));
    degopt=Degopt(graph);
    x=degopt.x;
    y=degopt.y;
    y .*= state.f(0)/eval_graph(graph,0);
    (graph,_)=graph_degopt(x,y);
    set_coeffs!(state.graph,get_coeffs(graph));
    err_after=eval_graph(graph,0)-state.f(0);
    println(" after $err_after");
    return state;
end

function runcommand(s::ShowInfo,state)

    graph=state.graph;
    try

        theta_guess=state.params[:rho];
        if (s.theta_guess != Nothing)
            theta_guess = s.theta_guess
        end

        (ff,theta)=compute_bwd_theta_exponential(graph,coefftype=BigFloat,
                                                tolerance=eps()/2,theta_init=big(theta_guess));

        @show theta
        @show ff(0)
        @show ff(1e-10)
        @show ff(theta_guess)
    catch e
        println("Unable to compute theta: $(s.theta_guess)");
    end

    return state;

end


function run_sequence(state,predefsims="");
    max_j=1000;
    statelist=Vector{Any}(undef,max_j);
    commandlist=Vector{String}(undef,max_j);
    statelist[1]=state;
    pre_commandlist=split(predefsims,"");

    err=showerr(state,output=false);
    println("Initialerr $(state.params[:graphname]):$err");
    j=1;
    command="";
    while (command != "q")
        state=statelist[j] # Base it on the old state
        print("Commandlist: "*String(join(commandlist[1:(j-1)]))* " ");

        x="";
        if (j <= size(pre_commandlist,1)) # It's stored in
            x=pre_commandlist[j];
            pre_commandlist[j] = "";
        end
        if x == ""
            x=String(read_one_key());
        end

        commandlist[j]=x;
        println("$x");
        if (x=="b") # Back
            j=j-2;
        elseif (x == "q")
            break
        else
            command = initcommand(x,state)
            newstate = runcommand(command,state)
            statelist[j+1] = newstate;
        end

        showerr(statelist[j+1]);

        j += 1;

    end

    state=statelist[j];
    return (state,statelist[1:j],join(commandlist[1:j]))

end


function showmetainfo(state)
    println("rho: $(state.params[:graphname]): ρ=$(state.params[:rho]) ");
end

function showerr(state;output=true,showmeta=false,n=state.params[:target_n])
    f=state.f;
    discr=get_disc_discr(state,n=n);


    graph=state.graph;

    err=norm((f.(discr)-eval_graph(graph,discr))./f.(discr),Inf)
    imagnorm=Float64(norm(imag.(get_coeffs(graph))));
    if (output)
        if (showmeta)
            showmetainfo(state)
        end

        if (imagnorm>0)
            imagstr="norm(imag(coeffs))= $imagnorm"
        else
            imagstr="";
        end
        println("Target error: $(Float64(err)) $imagstr");
    end

    return err
end


function read_one_key(; io = stdin)
    setraw!(raw) = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), io.handle, raw)
    setraw!(true)
    s=read(io, 1)
    setraw!(false)
    return s
end

## Some generation requires random "kick-starting". The following
## is a tabulation of the random numbers generated by
## Julia Version 1.7.0-DEV.1241 such that the generation
## can work also if the random number generators change in later versions.
mutable struct FakeRandomState
    k
end
function fake_randn(f::FakeRandomState,s)

    #Random.seed!(f.k);
    data=Dict(1=> [0.20154471885274503],
              2=> [0.20154471885274503;1.0340551364982935],
              3=> [0.20154471885274503;1.0340551364982935;-0.5035805877040641],
              28=> [0.359629402423705
 -0.06242791627632372
 -1.399567921508011
  0.4475847763985716
 -0.8809915030490731
 -0.2017627452986302
 -0.5526516215876736
 -0.35572577945591155
 -0.1607797466665553
  1.1475981747910613
  1.5035230388445053
  0.2812185565541045
 -1.6580378604283104
 -0.5724292571431209
  0.11417065387535505
 -0.1722032235811361
  0.10954150466894817
 -0.20957790746640031
  0.9978114970316212
 -0.33835189120333353
  0.32711196563255096
  1.3252610879493163
 -0.6792960939508351
 -0.3575966370924379
  2.1577149397345545
 -0.9170456413152324
 -1.021037073006077
                    -1.300754560220488],
              70=>
              [  0.359629402423705
 -0.06242791627632372
 -1.399567921508011
  0.4475847763985716
 -0.8809915030490731
 -0.2017627452986302
 -0.5526516215876736
 -0.35572577945591155
 -0.1607797466665553
  1.1475981747910613
  1.5035230388445053
  0.2812185565541045
 -1.6580378604283104
 -0.5724292571431209
  0.11417065387535505
 -0.1722032235811361
  0.10954150466894817
 -0.20957790746640031
  0.9978114970316212
 -0.33835189120333353
  0.32711196563255096
  1.3252610879493163
 -0.6792960939508351
 -0.3575966370924379
  1.8019224806064675
 -0.1963681193189925
 -0.3211451010184779
  1.2571922200069134
  1.7786053509127615
  0.14664320175311135
 -0.7581751321548498
  0.04258439124049732
  0.3352535968623635
 -0.5458913912138869
 -0.22338451264854198
 -0.34276021437488724
  0.3055002759810978
 -1.7936388425860077
 -1.3858822457781843
  0.19432728408125308
  0.9846027231525627
  0.15797132251418744
 -1.133614480683449
 -1.4271913217697576
  0.681897286305289
 -0.5715336031075239
 -0.4494921999102889
  1.85227626772481
 -0.7005686939516828
  1.2121103292882787
 -0.48418549713360953
  1.2973921022483264
  1.2098214008295365
 -0.9687971894606597
 -1.0648969106849027
  0.935342888589825
  0.8011405444777072
  0.6074278447599545
 -1.3830376414176226
  0.7946559924063759
 -0.9880332785507824
 -1.4752050514434305
 -0.5193534241313768
  2.048090887179751
  2.1577149397345545
 -0.9170456413152324
 -1.021037073006077
 -1.300754560220488
  0.10149275111695119
                 -0.11990415249480861],
              31 => [  0.359629402423705
 -0.06242791627632372
 -1.399567921508011
  0.4475847763985716
 -0.8809915030490731
 -0.2017627452986302
 -0.5526516215876736
 -0.35572577945591155
 -0.1607797466665553
  1.1475981747910613
  1.5035230388445053
  0.2812185565541045
 -1.6580378604283104
 -0.5724292571431209
  0.11417065387535505
 -0.1722032235811361
  0.10954150466894817
 -0.20957790746640031
  0.9978114970316212
 -0.33835189120333353
  0.32711196563255096
  1.3252610879493163
 -0.6792960939508351
 -0.3575966370924379
  2.1577149397345545
 -0.9170456413152324
 -1.021037073006077
 -1.300754560220488
  0.10149275111695119
 -0.11990415249480861
  0.7959522954870297]);

    #x=randn(s);
    #writedlm("/tmp/random_$(f.k)_$s.csv",x);
    x=data[s];
    f.k += s;
    return x;
end
