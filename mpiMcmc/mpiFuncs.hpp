#ifndef MPIFUNC_HPP
#define MPIFUNC_HPP

#include <array>
#include <fstream>

#include "Model.hpp"
#include "mpiMcmc.hpp"

void make_cholesky_decomp(struct ifmrMcmcControl &, Matrix<double, NPARAMS, nSave> &);
double logPostStep(Chain &, const Model &, std::array<double, N_WD_MASS1> &, Cluster &, double, std::array<double, 2> &, const std::vector<int>&);
int acceptClustMarg (double, double, std::array<double, 2> &);

void initChain (Chain &, const struct ifmrMcmcControl &, const Model &, std::array<double, 2> &, const std::vector<int>&);
void initIfmrMcmcControl (Cluster &, struct ifmrMcmcControl &, const Model &, Settings &);
void initMassGrids (std::array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &, std::array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &, std::array<double, N_WD_MASS1> &, const Chain &);

void propClustBigSteps (Cluster &, const struct ifmrMcmcControl &);
void propClustIndep (Cluster &, const struct ifmrMcmcControl &);
void propClustCorrelated (Cluster &, const struct ifmrMcmcControl &);

void printHeader (std::ofstream &, const std::array<double, NPARAMS> &);
void readCmdData (std::vector<Star> &, struct ifmrMcmcControl &ctrl, const Model &, std::vector<int>&);

#endif
