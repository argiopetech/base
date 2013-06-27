#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

#include "structures.hpp"

void evolve (struct cluster *the_cluster, std::vector<struct star> &stars, int index);
void calcAbsCoeffs (int filterSet);

#endif
