using GraphMatFun, LsqFit, LinearAlgebra
include("exp_reset_all.jl");
#include("simulationtools.jl");
include("simulationtools.jl");


f=exp
m=4
target_n=200 # Measure error wrt this discr
n=50;        # Optimize wrt this discr
rho=6.9e-1;  # SID rho m=4


its=5;
droptol0=1e-10;
opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                :linlsqr => :real_svd,:maxit => its)
base = init_state_mult(f,rho,n,m; eltype=Complex{BigFloat})
base.params[:n]=n;
base.params[:target_n]=target_n;
base.params[:opt_kwargs]=opt_kwargs;


# First the optimization free
sid_org=init_state_mult!(deepcopy(base),:sid,m,showmeta=true)

bbc_org=init_state_mult!(deepcopy(base),:bbc,m,showmeta=true)

sastre_org=init_state_mult!(deepcopy(base),:sastre,m,showmeta=true)

ps_org=init_state_mult!(deepcopy(base),:ps,m,showmeta=true)


mono_org=init_state_mult!(deepcopy(base),:mono,m,showmeta=true)



(mono,simlist,commandlist)=
     run_sequence(mono_org,"sGGGssddsddsddsddsggggggssssssGGGssssssssddsddsdsgggggggsdssssssssssssssssGGGsGGsGGsssssssdddsssssssssssssq");

(ps,simlist,commandlist)=
     run_sequence(ps_org,
     "sssssssssssddsssddssssdssGgggssdsssssssssssssDdgGGssssddssssssssssssssssddsssssssssssddssssssssssssddsssssssddsssssssssssssddgggssssssssssGssGGGssssssgggsdsssssssssssssssssssssssssddsssssssssssssssssssssssssssssssssddsssssssssssssssssssssssddssssddssssssddssddDssdggggsssssGssGGGGGGGsGsssssdsgggggggsddggsssssGGGGGGGsssssssssssssdsdsq")


(sastre,simlist,commandlist)=
     run_sequence(sastre_org,
                  "sssssddssssddsssssssssssdddssdddssssdddssdddsssNsssdddssssssssq");


(sid,simlist,commandlist)=
     run_sequence(sid_org,
     "sssssddssssddsssssssssssdddssdddssssdddssdddsssNsssdddsssssssssdddssssddsdsdssddsddsddsddsddsddsddsddsNssdsq")

(bbc,simlist,commandlist)=
     run_sequence(bbc_org,
     "ssssdddssssslddsllddssssssssssssddsssssssssddsssssssssdgggsssssssssssGGGsssssdddssssssssssq")


include("exp_print_all.jl");
println("Run include(exp_save_all.jl) if you want to save");
