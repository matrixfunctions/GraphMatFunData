using GraphMatFun
data_dir=joinpath("..","data","exp");
gen_code_dir=joinpath("..","data","generated","exp");
# Generated code
langlist=[LangJulia(),LangMatlab(), LangC_MKL(), LangC_OpenBLAS()];
extlist=[".jl",".m","_MKL.c","_OpenBLAS.c"];
for fname in readdir(data_dir);

    if (endswith(fname,".cgr"))
        # Only process cgr-files
        basename=split(fname,".")[1];
        println("Generating code for $fname");
        graph=import_compgraph(joinpath(data_dir,fname));
        compress_graph!(graph);

        if (eltype(graph) == BigFloat)
            graph=Compgraph(Float64,graph);
        end
        if (eltype(graph) == Complex{BigFloat})
            graph=Compgraph(ComplexF64,graph);
        end


        for (i,lang) in enumerate(langlist)
            ext=extlist[i];
            gen_code(joinpath(gen_code_dir,"$(basename)$(ext)"),
                     graph,
                     funname=basename,lang=lang);

        end

    end


end
