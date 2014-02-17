#ifndef MPIFUNC_HPP
#define MPIFUNC_HPP

#include <array>
#include <fstream>

#include "Model.hpp"
#include "mpiMcmc.hpp"

void make_cholesky_decomp(struct ifmrMcmcControl &, MVatrix<double, NPARAMS> &);

void printHeader (std::ofstream &, const std::array<double, NPARAMS> &);
std::vector<StellarSystem> readCmdData (struct ifmrMcmcControl &ctrl, const Model &, std::vector<int>&, std::array<double, FILTS>&, std::array<double, FILTS>&, const Settings&);

#endif
