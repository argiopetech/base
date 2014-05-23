#ifndef MPIFUNC_HPP
#define MPIFUNC_HPP

#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "Model.hpp"
#include "mpiMcmc.hpp"

void printHeader (std::ofstream &, const std::array<double, NPARAMS> &);

#endif
