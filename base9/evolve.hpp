#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

#include "structures.hpp"

void evolve (Cluster *the_cluster, std::vector<Star> &stars, int index);
void calcAbsCoeffs (int filterSet);

#endif
