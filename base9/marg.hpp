#ifndef MARG_HPP
#define MARG_HPP

#include <vector>

#include "Star.hpp"
#include "Cluster.hpp"

#include "Model.hpp"
#include "Utility.hpp"

std::vector<double> margEvolveWithBinary (const Cluster &, std::vector<StellarSystem>&, const Model&, const Isochrone&, base::utility::ThreadPool&, const bool);
std::vector<double> margEvolveNoBinaries(const Cluster&, const Model&, const Isochrone&, base::utility::ThreadPool&, const double* const, const double* const, const double* const, const size_t, const size_t, const size_t, const bool);

#endif
