#ifndef EVOLVE_H
#define EVOLVE_H

#include <array>
#include <vector>

#include "Cluster.hpp"
#include "Star.hpp"

#include "Model.hpp"
#include "structures.hpp"

std::array<double, FILTS> evolve (const Cluster&_cluster, Model const&, const std::vector<int>&, StellarSystem&);

#endif
