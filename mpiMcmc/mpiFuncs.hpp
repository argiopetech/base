#ifndef MPIFUNC_HPP
#define MPIFUNC_HPP

#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "Model.hpp"
#include "mpiMcmc.hpp"

void printHeader (std::ofstream &, const std::array<double, NPARAMS> &);
std::pair<std::vector<std::string>, std::vector<StellarSystem>> readCmdData (struct ifmrMcmcControl &ctrl, std::vector<double>&, std::vector<double>&, const Settings&);

#endif
