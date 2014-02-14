#ifndef EVOLVE_H
#define EVOLVE_H

#include <array>
#include <vector>

#include "Cluster.hpp"
#include "Star.hpp"

#include "Model.hpp"
#include "structures.hpp"

std::array<double, FILTS> evolve (const Cluster &the_cluster, Model const &, const std::vector<int>&, Star &star, std::array<double, 2>&);

#endif
