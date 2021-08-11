# GraphMatFunData

This repository contains data files for the julia package [`GraphMatFun.jl`](https://github.com/matrixfunctions/GraphMatFun.jl) in particular the data files for the experiments in the manuscript [Computation graph for matrix fucntions](https://arxiv.org/abs/2107.12198).

## Examples:


### Julia

You can use the generated files to efficiently compute the matrix exponential.

```julia
julia> using GraphMatFun, BenchmarkTools, LinearAlgebra;
julia> # Download from https://github.com/matrixfunctions/GraphMatFunData/tree/main/data/generated/exp
julia> include("exp_sid_m5_opt_rho1_9.jl");
julia> A=randn(1000,1000); A=1.8*A/norm(A);
julia> E1=@btime exp($A);
  1.111 s (23 allocations: 76.31 MiB)
julia> E2=@btime exp_sid_m5_opt_rho1_9($A);
  332.667 ms (19 allocations: 68.67 MiB)
julia> norm(E1-E2)/norm(E2)
3.1592474046054336e-16
```

### Matlab

The following illustrates that the fixed topology optimization improves the error.

```matlab
>> % Download exp_sid_m5.m and exp_sid_m5_opt_rho1_9.m from https://github.com/matrixfunctions/GraphMatFunData/tree/main/data/generated/exp
>> xv=0:0.01:3;
>> err1=zeros(size(xv)); for j=1:length(xv); err1(j)=exp_sid_m5(xv(j))-exp(xv(j)); end
>> err2=zeros(size(xv)); for j=1:length(xv); err2(j)=exp_sid_m5_opt_rho1_9(xv(j))-exp(xv(j)); end
>> semilogy(xv,abs(err1)./exp(xv))
>> clf
>> semilogy(xv,abs(err1)./exp(xv))
>> hold on
>> semilogy(xv,abs(err2)./exp(xv))
>> pwd
```

Produces:

 ![Screenshot_20210806_102931](https://user-images.githubusercontent.com/11163595/128481669-dfbb7a85-7e20-4d9c-8dcb-937eba9954d4.png)


## Directory structure:
* `data/exp`: The graphs (in CGR-format) for the exponential
* `data/sqrt_plus_one`: The graphs (in CGR-format) for the square root
* `data/generated`: The code generated (in matlab, julia and C-format) for all the filies in `data/exp` and `data/sqrt_plus_one`
* `src/opt`: The design by optimization of the exponential and square root. Can generate files in `data/`
* `src/cputimg`: The cpu-timing for the expontial. Reads the data files from `data/`
* `src/matlab`: Loading and evaluation of degree optimal polynomials in matlab
