using MD5
function restricted_md5(fname)
    if (!isfile(fname))
        return randn(10);
    end
    f=open(fname,"r");
    lines=readlines(fname);
    for (i,line) in enumerate(lines)
        if occursin("Created: ", line)
            lines[i]=""
        end
    end
    md5(join(lines,"\n"));
end

include("exp_set_list_all.jl");


println("Saving simulation");

rho=mono.params[:rho];
m=count(values(state_list[1].graph.operations) .== :mult)

println("Domain: $(rho), m=$m",);
for (i,state)=enumerate(state_list)
    if (!isnothing(state))
        name=state.params[:graphname]
        @show name
        if contains(name,"_opt")
            graph=state.graph;
            rho_str=replace(string(rho),"."=>"_");
            fname=joinpath("..","..","data","exp","exp_$(name)_m$(m)_$(rho_str).cgr");
            err=Float64(showerr(state,output=false,n=state.params[:target_n]*2))

            isreal=false;
            if (norm(imag.(get_coeffs(graph)))==0)
                isreal=true;
                g=Compgraph(real(eltype(graph)),graph);
            end
            println("Saving $name $fname isreal=$isreal");
            optimized=""
            if (contains(name,"+GN"))
                optimized=" optimized starting from "*replace(mono.params[:graphname], "+GN" => "");
            end


            tmpname = "/tmp/tmpfile.jl"
            export_compgraph(g, tmpname;
                             fun="exp", dom="$(rho)D", err="$err",
                             order=get_topo_order(g)[1],
                             descr="Matrix exponential with degree optimal form $optimized",
                             user="Elias Jarlebring");

            # Only save it if the file is actually different (beside time-stamp)
            if (restricted_md5(tmpname) != restricted_md5(fname))
                cp(tmpname,fname);
            else
                println("Skip saving $i (since file unchanged)");
            end

        else
            println("Skip saving $i (since it is the non-optimized)");
        end

    else
        println("Skip saving $i");
    end
end
