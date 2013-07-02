#ifndef MPIFUNC_HPP
#define MPIFUNC_HPP

#include <array>

#include "Model.hpp"
#include "mpiMcmc.hpp"

void make_cholesky_decomp(struct ifmrMcmcControl &ctrl, Matrix<double, NPARAMS, nSave> &params);
double logPostStep(Chain &mc, const Model &evoModel, std::array<double, N_WD_MASS1> &wdMass1Grid, Cluster &propClust, double fsLike, std::array<double, 2> &ltau);
int acceptClustMarg (double logPostCurr, double logPostProp, std::array<double, 2> &ltau);
void initChain (Chain &mc, const struct ifmrMcmcControl &ctrl, const Model &evoModels, std::array<double, 2> &ltau);

#endif
