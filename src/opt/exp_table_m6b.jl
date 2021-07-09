using GraphMatFun, LsqFit, LinearAlgebra, GenericSVD,Random
include("exp_reset_all.jl");
#include("simulationtools.jl");
include("simulationtools.jl");


f=exp
m=6
target_n=200 # Measure error wrt this discr
n=50;        # Optimize wrt this discr
rho=2.7;  # Increased SID rho m=5


its=5;
droptol0=1e-10;
opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                :linlsqr => :real_svd,:maxit => its)
base = init_state_mult(f,rho,n,m; eltype=Complex{BigFloat})
base.params[:n]=n;
base.params[:target_n]=target_n;
base.params[:opt_kwargs]=opt_kwargs;
base.params[:kickit_mode]=4;


oldrho="2_22";
datadir=joinpath("..","..","data","exp");
# First the optimization free
filename=joinpath(datadir,"exp_sid_m$(m).cgr");
sid_org=init_state_file!(deepcopy(base),:sid,filename,showmeta=true)

#filename=joinpath("exp_bbc_m$(m)_rho$oldrho.cgr");
#bbc_org=init_state_file!(deepcopy(base),:bbc,filename,showmeta=true)

filename=joinpath(datadir,"exp_sastre_m$(m).cgr");
sastre_org=init_state_file!(deepcopy(base),:sastre,filename,showmeta=true)


#sastre_org.cref=get_degopt_crefs(sastre_org.graph);


filename=joinpath(datadir,"exp_ps_m$(m)_rho$oldrho.cgr");
ps_org=init_state_file!(deepcopy(base),:ps,filename,showmeta=true)


filename=joinpath(datadir,"exp_mono_m$(m)_rho$oldrho.cgr");
mono_org=init_state_file!(deepcopy(base),:mono,filename,showmeta=true)


filename=joinpath(datadir,"exp_mono_m6_opt_rho2_7_warmstart.cgr");
#filename="/home/jarl/jobb/src/matfun/simulations/newgraphs/exp_m$(m)_mono+GN_$oldrho.cgr";
mono0=init_state_file!(deepcopy(base),:mono,filename,showmeta=true)
(mono,simlist,commandlist)=
     run_sequence(mono0,"dddddddddddddddddddddsssGGGssdddssssdddddsssssssdddsdddsssszq");
#"dddddddddddddddddddddsssGGGSSSddddSkSSdSKKKKKKsSdddgggggggggSSSSSSSSsssssssssssSSDsSSSSSSSSSSDSSSSGGGGGGGGGSSSdSdkDSSztkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkSSq"
filename=joinpath(datadir,"exp_ps_m$(m)_opt_rho$oldrho.cgr");
ps0=init_state_file!(deepcopy(base),:ps,filename,showmeta=true)
ps0.params[:kickit_mode]=2;
(ps,simlist,commandlist)=
     run_sequence(ps0,
     "sdddddddddSdddSdSddddSSNSkSSSGGGSddSddSSddSSzddSddkSkkkkkkkkkkkSdSgggggggSSSSSSSSssGGGGGGGSSSzSddssssddsddszk");


filename=joinpath(datadir,"exp_sastre_m$(m)_opt_rho$oldrho.cgr");
sastre0=init_state_file!(deepcopy(base),:sastre,filename,showmeta=true)
(sastre,simlist,commandlist)=
     run_sequence(sastre0,
                  "ddddddddddsdddddddddddssGGGssddddsdddddsssdddssszq");

filename=joinpath(datadir,"exp_sid_m$(m)_opt_rho$oldrho.cgr");
sid0=init_state_file!(deepcopy(base),:sid,filename,showmeta=true)
(sid,simlist,commandlist)=
     run_sequence(sid0,
     "ddddddddddsdddddddddddssGGGssddddsdddddsssdddsssdddsssddddsssssdsddszq")

#filename=joinpath("exp_bbc+GN_m$(m)_rho$oldrho.cgr");
#bbc0=init_state_file!(deepcopy(base),:bbc,filename,showmeta=true)
#(bbc,simlist,commandlist)=
#     run_sequence(bbc0,
#     "sdddddsdddddssssssddddsssssssddddsssssssssssssdddssssssssssddddssssssssdsggggggsdssssssssssssssssGGGGGsssGGGsssq")
#

include("exp_print_all.jl");
println("Run include(exp_save_all.jl) if you want to save");
