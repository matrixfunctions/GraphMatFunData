Matrix exponential approximations:

## Basic
* `exp_mono_mX.cgr`: monomial approximation with X multiplications (available from `GraphMatFun.graph_monomial`) with Taylor coefficients
* `exp_mono_mX_rhoZ_Z.cgr`: monomial approximation with X multiplications (available from `GraphMatFun.graph_monomial`) with coefficients determined from the linead least squares fit on the disk with radius Z.Z.
* `exp_ps_mX.cgr`: Paterson-Stockmayer polynomial evaluation with X multiplications (available from `GraphMatFun.graph_ps`) with Taylor coefficients
* `exp_mono_mX_rhoZ_Z.cgr`: Paterson-Stockmayer polynomial evaluation with X multiplications (available from `GraphMatFun.graph_ps`) with coefficients determined from the linead least squares fit on the disk with radius Z.Z.
* `exp_sastre_mX.cgr`: Sastre paper paper: with X multiplications (implemented in `GraphMatFun.graph_sastre_exp`)
* `exp_sid_mX.cgr`: SID paper paper: with X multiplications (implemented in `GraphMatFun.graph_sid`)
* `exp_bbc_mX.cgr`: BBC paper: with X multiplications (implemented in `GraphMatFun.graph_bbc_exp`)
* `exp_bbcs_mX.cgr`: BBCS paper: with X multiplications (implemented in `GraphMatFun.graph_bbcs_exp`)
* `exp_native_jl_rhoZ_Z.cgr`: The method implemented in julias function `exp(A)` suitable for matrices with norm less than Z.Z

## Optimized 
* `exp_mono_mX_opt_rhoZ_Z.cgr`: the optimization of `exp_mono_mX` in degopt format on the disk with radius rho=Z.Z
* `exp_ps_mX_opt_rhoZ_Z.cgr`: the optimization of `exp_ps_mX` in degopt format on the disk with radius rho=Z.Z
* `exp_sastre_mX_opt_rhoZ_Z.cgr`: the optimization of `exp_sastre_mX` in degopt format on the disk with radius rho=Z.Z
* `exp_sid_mX.cgr`: the optimization of `exp_sid_mX` in degopt format on the disk with radius rho=Z.Z
* `exp_bbc_mX_opt_rhoZ_Z.cgr`: the optimization of `exp_bbc_mX` in degopt format on the disk with radius rho=Z.Z
