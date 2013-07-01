#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

#include "Model.hpp"
#include "structures.hpp"

void evolve (Cluster *the_cluster, Model &, Star &star);
void calcAbsCoeffs (int filterSet);

#endif
