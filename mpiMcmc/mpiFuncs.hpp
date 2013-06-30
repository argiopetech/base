#include <array>

#include "Model.hpp"
#include "mpiMcmc.hpp"

void make_cholesky_decomp(struct ifmrMcmcControl &ctrl, Matrix<double, NPARAMS, nSave> &params);
double logPostStep(Chain &mc, Model &evoModel, std::array<double, N_WD_MASS1> &wdMass1Grid, Cluster &propClust, double fsLike);
int acceptClustMarg (double logPostCurr, double logPostProp);
