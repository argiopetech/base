#include <array>

#include "mpiMcmc.hpp"

void make_cholesky_decomp(struct ifmrMcmcControl &ctrl, double **params);
double logPostStep(struct chain &mc, std::array<double, N_WD_MASS1> &wdMass1Grid, struct cluster &propClust, double fsLike);
int acceptClustMarg (double logPostCurr, double logPostProp);
