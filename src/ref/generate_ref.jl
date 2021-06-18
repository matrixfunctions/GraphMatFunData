# export JULIA_LOAD_PATH=$JULIA_LOAD_PATH:~/jobb/src/matfun
using GraphMatFun
mv=1:10
data_dir=joinpath("..","..","data","exp");
references=Dict{Symbol,String}();
references[:sid]="Boosting the computation of the matrix exponential, J. Sastre, J. Ibáñez, E. Defez, Appl. Math.  Computation, 340, 2019, 206-220";
references[:sastre]="Efficient evaluation of matrix polynomials, J. Sastre. Linear Algebra and its Applications ,Volume 539, 2018, Pages 229-250, https://doi.org/10.1016/j.laa.2017.11.010";
references[:bbc]="Computing the matrix exponential with an optimized Taylor polynomial approximation, P. Bader, S.  Blanes, and F. Casas, Mathematics, 7(12), 2019";
references[:ps]="On the number of nonscalar multiplications necessary to evaluate polynomials, M. Paterson, L.   Stockmeyer, SIAM J. Comput., 2(1), 1973";

methods=keys(references);


for m in mv
    for method in methods
        reftext=references[method]
        try
            extra_text="";
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
            end
            compress_graph!(graph);
            descr="Reference implementation of $reftext with m=$m multiplications$extra_text"
            export_compgraph(graph,joinpath(data_dir,"$(method)_m$m.cgr"),
                             descr=descr,user="Elias Jarlebring",
                             genby="GraphMatFunData/src/ref/generate_ref.jl");
        catch (e)
            println("$method: Skipping m=$m multiplications");
        end
    end

end
