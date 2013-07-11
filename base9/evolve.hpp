#ifndef EVOLVE_H
#define EVOLVE_H

#include <array>
#include <vector>

#include "Model.hpp"
#include "structures.hpp"

void evolve (const Cluster &the_cluster, Model const &, std::array<double, FILTS>&, const std::vector<int>&, Star &star, std::array<double, 2>&);

#endif
