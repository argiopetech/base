#include <array>

#include "mpiMcmc.hpp"

void make_cholesky_decomp(struct ifmrMcmcControl &ctrl, Matrix<double, NPARAMS, nSave> &params);
double logPostStep(struct chain &mc, std::array<double, N_WD_MASS1> &wdMass1Grid, Cluster &propClust, double fsLike);
int acceptClustMarg (double logPostCurr, double logPostProp);
