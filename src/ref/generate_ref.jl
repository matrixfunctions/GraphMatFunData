# export JULIA_LOAD_PATH=$JULIA_LOAD_PATH:~/jobb/src/matfun
using GraphMatFun
mv=1:10
data_dir=joinpath("..","..","data","exp");
references=Dict([:sid => "Boosting the computation of the matrix exponential, J. Sastre, J. Ibáñez, E. Defez, Appl. Math.  Computation, 340, 2019, 206-220",
                :sastre => "Efficient evaluation of matrix polynomials, J. Sastre. Linear Algebra and its Applications ,Volume 539, 2018, Pages 229-250, https://doi.org/10.1016/j.laa.2017.11.010"]);
methods=[:sid,:sastre];

for m in mv
    for method in methods
        reftext=references[method]
        try
            if (method == :sid)
                (graph,_)=graph_sid_exp(m);
            elseif (method == :sastre)
                (graph,_)=graph_sastre_basic_exp(m);
            end
            compress_graph!(graph);
            descr="Reference implementation of $reftext with m=$m multiplications"
            export_compgraph(graph,joinpath(data_dir,"$(method)_m$m.cgr"),
                             descr=descr,user="Elias Jarlebring",
                             genby="GraphMatFunData/src/ref/generate_ref.jl");
        catch (e)
            println("$method: Skipping m=$m multiplications");
        end
    end

end
