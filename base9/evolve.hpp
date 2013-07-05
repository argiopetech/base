#ifndef EVOLVE_H
#define EVOLVE_H

#include <array>
#include <vector>

#include "Model.hpp"
#include "structures.hpp"

void evolve (Cluster &the_cluster, Model const &, Star &star, std::array<double, 2>&);

#endif
