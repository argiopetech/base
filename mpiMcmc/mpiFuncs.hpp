#ifndef MPIFUNC_HPP
#define MPIFUNC_HPP

#include <array>
#include <fstream>

#include "Model.hpp"
#include "mpiMcmc.hpp"

void make_cholesky_decomp(struct ifmrMcmcControl &, MVatrix<double, NPARAMS> &);

void initChain (Chain &, const Model &, std::array<double, 2> &, const std::vector<int>&);

void printHeader (std::ofstream &, const std::array<double, NPARAMS> &);
void readCmdData (std::vector<Star> &, struct ifmrMcmcControl &ctrl, const Model &, std::vector<int>&, std::array<double, FILTS>&, std::array<double, FILTS>&, const Settings&);

#endif
