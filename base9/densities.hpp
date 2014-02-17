/* densities.h */

#ifndef DENSITIES_H
#define DENSITIES_H

#include <array>
#include <vector>

#include "Cluster.hpp"
#include "constants.hpp"
#include "Model.hpp"
#include "Star.hpp"

double logPriorClust (const Cluster&, const Model&);
double logPost1Star (const StellarSystem&, const Cluster&, const Model&, const std::vector<int>&);
double logTDens (double, double, double, double);
#endif
