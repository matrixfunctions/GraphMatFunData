# export JULIA_LOAD_PATH=$JULIA_LOAD_PATH:~/jobb/src/matfun
using GraphMatFun
mv=1:11
data_dir=joinpath("..","..","data","exp");
references=Dict{Symbol,String}();
references[:sid]="Boosting the computation of the matrix exponential, J. Sastre, J. Ibáñez, E. Defez, Appl. Math.  Computation, 340, 2019, 206-220";
references[:sastre]="Efficient evaluation of matrix polynomials, J. Sastre. Linear Algebra and its Applications ,Volume 539, 2018, Pages 229-250, https://doi.org/10.1016/j.laa.2017.11.010";
references[:bbc]="Computing the matrix exponential with an optimized Taylor polynomial approximation, P. Bader, S.  Blanes, and F. Casas, Mathematics, 7(12), 2019";
references[:ps]="On the number of nonscalar multiplications necessary to evaluate polynomials, M. Paterson, L.   Stockmeyer, SIAM J. Comput., 2(1), 1973";
references[:mono]="monomial Taylor series evaluation";
references[:exp_native_jl]="converted Julia 1.7 implementation of Scaling and squaring https://github.com/JuliaLang/julia/blob/697e782ab86bfcdd7fd15550241fe162c51d9f98/stdlib/LinearAlgebra/src/dense.jl#L554 which is based on N. J. Higham. The Scaling and Squaring Method for the Matrix Exponential Revisited. SIAM J. Matrix Anal. Appl., 2005 26:4, 1179-1193"
methods=keys(references);


for m in mv
    for method in methods
        reftext=references[method]
        multstr="with m=$m multiplications";
        fname="$(method)_m$m.cgr";
        try
            extra_text="";
            dom=""; # Domain
            if (method == :sid)
                (graph,_)=graph_sid_exp(m);
            elseif (method == :sastre)
                (graph,_)=graph_sastre_basic_exp(m);
            elseif (method == :bbc)
                if (m==6)
                    println("Skipping $method m=6");
                    continue;
                end
                (graph,_)=graph_bbc_basic_exp(m);

            elseif (method == :ps)
                # Automatically determine degree
                deg=0;
                old_graph=NaN; graph=NaN;
                for deg0=1:100
                    taylor=1 ./ factorial.(convert.(BigFloat,0:deg0));
                    old_graph=graph;
                    (graph,_)=graph_ps(taylor);
                    if (count(values(graph.operations) .== :mult) > m)
                        deg=deg0-1;
                        break;
                    end
                end
                graph=old_graph;
                extra_text=" degree=$deg";
            elseif (method == :mono)
                # Automatically determine degree
                deg=0;
                old_graph=NaN; graph=NaN;
                for deg0=1:100
                    taylor=1 ./ factorial.(convert.(BigFloat,0:deg0));
                    old_graph=graph;
                    (graph,_)=graph_monomial_degopt(taylor);
                    if (count(values(graph.operations) .== :mult) > m)
                        deg=deg0-1;
                        break;
                    end
                end
                graph=old_graph;
                extra_text=" degree=$deg";
            elseif method == :exp_native_jl
                normlist=[0.015;0.25;0.95;2.1;5.4];

                for j=1:m-size(normlist,1)
                    push!(normlist,normlist[end]*2); # square
                end



                n=normlist[m] # Interpret m as the m'th rho interval
                (graph,_)=graph_exp_native_jl(reshape([n-0.0001],1,1));
                rho=n;
                dom="$(rho)D";
                mults=count(values(graph.operations) .== :mult)
                multstr="with $mults multiplications and one inverse multiplication";
                rhostr=replace(string(rho),"." => "_");
                fname="exp_native_jl_rho$(rhostr).cgr";
            end
            compress_graph!(graph);
            descr="Reference implementation of $reftext $multstr$extra_text"
            export_compgraph(graph,joinpath(data_dir,fname),
                             descr=descr,user="Elias Jarlebring",
                             dom=dom,
                             genby="GraphMatFunData/src/ref/generate_ref.jl");
        catch (e)
            if (typeof(e) != ErrorException)
                rethrow(e);
            end

            println("$method: Skipping m=$m multiplications");
        end
    end

end

gen_code_dir=joinpath("..","..","data","generated","exp");
# Generated code
langlist=[LangJulia(),LangMatlab(), LangC_MKL(), LangC_OpenBLAS()];
extlist=[".jl",".m","_MKL.c","_OpenBLAS.c"];
for fname in readdir(data_dir);

    if (endswith(fname,".cgr"))
        # Only process cgr-files
        basename=split(fname,".")[1];
        @show basename
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
