using GraphMatFun, LsqFit, LinearAlgebra, GenericSVD,Random
include("exp_reset_all.jl");
#include("simulationtools.jl");
include("simulationtools.jl");


f=exp
m=8
target_n=400 # Measure error wrt this discr
n=100;        # Optimize wrt this discr
rho=13.5;  # Increased SID rho m=7


its=5;
droptol0=1e-20;
opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                :linlsqr => :real_svd,:maxit => its,
                :errtype => :relerr)
base = init_state_mult(f,rho,n,m; eltype=Complex{BigFloat})
base.params[:n]=n;
base.params[:target_n]=target_n;
base.params[:opt_kwargs]=opt_kwargs;
base.params[:kickit_mode]=1;
base.params[:cref_count]=0;

oldrho="6_4";
oldrho2="3_59";
## First the optimization free
#filename="simulations/newgraphs/exp_m$(m)_sid_$oldrho.cgr";
#sid_org=init_state_file!(deepcopy(base),:sid,filename,showmeta=true)
#
##filename="simulations/newgraphs/exp_m$(m)_bbc_$oldrho.cgr";
##bbc_org=init_state_file!(deepcopy(base),:bbc,filename,showmeta=true)
#
#filename="simulations/newgraphs/exp_m$(m)_sastre_$oldrho.cgr";
#sastre_org=init_state_file!(deepcopy(base),:sastre,filename,showmeta=true)
#
#
##sastre_org.cref=get_degopt_crefs(sastre_org.graph);
#
#
#filename="simulations/newgraphs/exp_m$(m)_ps_$oldrho.cgr";
#ps_org=init_state_file!(deepcopy(base),:ps,filename,showmeta=true)
#
#
#filename="simulations/newgraphs/exp_m$(m)_mono_$oldrho.cgr";
#mono_org=init_state_file!(deepcopy(base),:mono,filename,showmeta=true)
#


#filename="simulations/newgraphs/exp_m$(m-1)_mono+GN_$oldrho.cgr";
##filename="simulations/graphs/exp_m8_mono_taylor_13_5.cgr";
#mono0=init_state_file!(deepcopy(base),:mono,filename,showmeta=true,scale_and_square=true)
#mono0.cref=mono0.cref[8:70]
#(mono,simlist,commandlist)=
#     run_sequence(mono0,"ksssssdddsssssddddssssssddddddsssssssdddsssssssssssdddssssGGsssddddggggggssssssssGGGsssGGGssssdddddggggggggsssssssGGGsGGGGGssssdddddggggggggssssssssGGGsGGGGGsssssddddgggggggsssGGGGGGGssssssdddggggggggggssssssggssssssGsGsGsGsGGsGGsGGGsGGsss");
#

datadir=joinpath("..","..","data","exp");
#filename="simulations/graphs/exp_m8_mono_taylor_13_5.cgr";
filename=joinpath(datadir,"exp_mono_m$(m-1)_opt_rho$(oldrho).cgr");
mono0=init_state_file!(deepcopy(base),:mono,filename,showmeta=true,
                       scale_and_square=true)

degopt_crefs=get_degopt_crefs(mono0.graph);
#mono0.cref=mono0.cref[8:end-8]
mono0.cref=[#get_degopt_crefs(mono0.graph)[1][end-1][1];
            #get_degopt_crefs(mono0.graph)[1][end-1][2];
            degopt_crefs[1][end][1];
            degopt_crefs[1][end][2];
            degopt_crefs[2] ];
(mono,simlist,commandlist)=
run_sequence(deepcopy(mono0),"sSSNSkSSGGGSS22SdddddSSSSSzS0pSZpSSSSzSzq");




filename=joinpath(datadir,"exp_ps_m$(m-1)_opt_rho$(oldrho).cgr");
ps0=init_state_file!(deepcopy(base),:ps,filename,showmeta=true,
                         scale_and_square=true)
ps0.cref=[degopt_crefs[1][end][1];
            degopt_crefs[1][end][2];
            degopt_crefs[2]];

(ps,simlist,commandlist)=
     run_sequence(ps0,
                  "SSNSSSSGGGSSSSSDDDSSSSZSSSZSSSddSddSzZSSSSddZSSSSZpSSSzq");

filename=joinpath(datadir,"exp_sastre_m$(m-1)_opt_rho$(oldrho).cgr");
sastre0=init_state_file!(deepcopy(base),:sastre,filename,showmeta=true,
                         scale_and_square=true)
sastre0.cref=[degopt_crefs[1][end][1];
            degopt_crefs[1][end][2];
            degopt_crefs[2]];

(sastre,simlist,commandlist)=
     run_sequence(sastre0,
                  "SSSSNSSSSkSSGGGSSddddSdSdSzSkZSSSSSSSZsSSzZSSSSpZpSSSSddSddddSkSSS2zq");



#filename=joinpath(datadir,"exp_sid_m$(m-1)_opt_rho$(oldrho).cgr");
#sid0=init_state_file!(deepcopy(base),:sid,filename,showmeta=true,
#                         scale_and_square=true)
#sid0.cref=[degopt_crefs[1][end][1];
#            degopt_crefs[1][end][2];
#            degopt_crefs[2]];
#
#(sid,simlist,commandlist)=
#     run_sequence(sid0,
#                  "SSSSNSSSSkSSGGGSSddddSdSdSzSkZSSSSSSSZsSSzZSSSSpZpSSSSddSddddSkSSS2zq");
#
#asd
#sSSNSkSSGGGSS22Sq");
##sSkSSSddddSSSSSSS22GGG22ddd2SS22n222tk22222222223tdddtSS333k333333333dddd3N333k3333333k3333333k3333333dddds333333333333333333333sSSStztkk3333N3nn3333n33NN333n33tz33NNN3333N3tS333ts333t33Skkk3333333tTz");
#
#graph1=mono1.graph;
##degopt=Degopt(mono1.graph);
###degopt.x[1][1] += sqrt(eps());;
##degopt.x[1][1][1] += 0.000001;
##normalize!(degopt);
##degopt.x[1][1][1] -= 0.000001;
##graph1=graph_degopt(degopt)[1];
#export_compgraph(graph1,"/tmp/mono1_m8_tmp.cgr");
#
#
#mono0.cref=[#get_degopt_crefs(mono0.graph)[1][end-1][1];
#            #get_degopt_crefs(mono0.graph)[1][end-1][2];
#            get_degopt_crefs(mono0.graph)[1][end][1];
#            get_degopt_crefs(mono0.graph)[1][end][2];
#            get_degopt_crefs(mono0.graph)[2] ];
#filename_tmp="/tmp/mono1_m8_tmp.cgr";
#mono2=init_state_file!(deepcopy(base),:mono,filename_tmp,showmeta=true);
#setdiff!(mono2.cref,mono0.cref);
#(mono,simlist,commandlist)=
#run_sequence(deepcopy(mono2),"dddddddddddddddddd22GGG22n2S");
#
#
#
#

#(mono,simlist,commandlist)=
#run_sequence(deepcopy(mono0),"ssssGGGssgddssssssddddggggggggggsssGGGsss");


#             "kssssGGGssgddssssssddddggggggggggsssGGGssssGGGssGGGGsssssksssdddsssssddssssssdddsdgggggggsssssGGGGsGGGssssssNssGssddssssddssssdggggggggsssssGGGGssssGGGGsssddsdsdgggggggggggssssGGGGGGGGGGsGssssggggggggdsdsdskssGGGGGGGGssNsszq");
#sGsGsGGGssssssssddsssGGGGGGsssksssssNsq
#             "sssGGGsssdddsddsssssdddsssdddsssssq");

#filename="simulations/newgraphs/exp_m$(m-1)_mono+GN_$(oldrho).cgr";
#mono0=init_state_file!(deepcopy(base),:mono,filename,showmeta=true,scale_and_square=true)
#(mono,simlist,commandlist)=
#     run_sequence(mono0,"ddddddddddddddddddddddddddsssGGGsssdddsddsssssdddsssdddsssssq");
#
#filename="simulations/newgraphs/exp_m$(m)_ps+GN_$oldrho.cgr";
#ps0=init_state_file!(deepcopy(base),:ps,filename,showmeta=true)
#(ps,simlist,commandlist)=
#     run_sequence(ps0,
#     "sq");
#
#filename="simulations/newgraphs/exp_m$(m-1)_sastre+GN_$oldrho.cgr";
#sastre0=init_state_file!(deepcopy(base),:sastre,filename,showmeta=true,scale_and_square=true)
#(sastre,simlist,commandlist)=
#     run_sequence(sastre0,
#                  "kssGGGssddddsdddddsdsssssdddssssddssssddsgggsssssGGGssq");
#
#filename="simulations/newgraphs/exp_m$(m)_sid+GN_$oldrho.cgr";
#sid0=init_state_file!(deepcopy(base),:sid,filename,showmeta=true)
#(sid,simlist,commandlist)=
#     run_sequence(sid_org,
#     "nssNdddddddddddddddddddssGGGsssddsddsddsddsd")
#
##filename="simulations/newgraphs/exp_m$(m)_bbc+GN_$oldrho.cgr";
##bbc0=init_state_file!(deepcopy(base),:bbc,filename,showmeta=true)
##(bbc,simlist,commandlist)=
##     run_sequence(bbc0,
##     "sdddddsdddddssssssddddsssssssddddsssssssssssssdddssssssssssddddssssssssdsggggggsdssssssssssssssssGGGGGsssGGGsssq")
##
#
include("exp_print_all.jl");
println("Run include(exp_save_all.jl) if you want to save");
