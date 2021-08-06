# GraphMatFunData
Data files and scripts for GraphMatFun.jl 

This file contains data files for the julia package [`GraphMatFun.jl`](https://github.com/matrixfunctions/GraphMatFun.jl) in particular the data files for the experiments in the manuscript [Computation graph for matrix fucntions](https://arxiv.org/abs/2107.12198)

Directory structure:
* `data/exp`: The graphs (in CGR-format) for the exponential
* `data/sqrt_plus_one`: The graphs (in CGR-format) for the square root
* `data/generated`: The code generated (in matlab, julia and C-format) for all the filies in `data/exp` and `data/sqrt_plus_one`
* `src/opt`: The design by optimization of the exponential and square root. Can generate files in `data/`
* `src/cputimg`: The cpu-timing for the expontial. Reads the data files from `data/`
