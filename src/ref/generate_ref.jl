# export JULIA_LOAD_PATH=$JULIA_LOAD_PATH:~/jobb/src/matfun
using GraphMatFun
mv=[1;2;3;4;5;6;7]
data_dir=joinpath("..","..","data","exp");
references=Dict([:sid => "Boosting the computation of the matrix exponential, J. Sastre, J. Ibáñez, E. Defez, Appl. Math.  Computation, 340, 2019, 206-220"]);
methods=[:sid];

for m in mv
    for method in methods
        try
            (graph,_)=graph_sid_exp(m);
            compress_graph!(graph);
            reftext=references[method]
            descr="Reference implementation of $reftext with m=$m multiplications"
            export_compgraph(graph,joinpath(data_dir,"$(method)_m$m.cgr"),
                             descr=descr,user="Elias Jarlebring",
                             genby="GraphMatFunData/src/ref/generate_ref.jl");
        catch (e)
            println("$method: Skipping m=$m multiplications");
        end
    end

end
