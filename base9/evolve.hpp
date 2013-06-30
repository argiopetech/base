#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

#include "Model.hpp"
#include "structures.hpp"

void evolve (Cluster *the_cluster, Model &, std::vector<Star> &stars, int index);
void calcAbsCoeffs (int filterSet);

#endif
