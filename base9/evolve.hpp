#ifndef( EVOLVE_H )
#define EVOLVE_H

#include "structures.hpp"

void evolve (struct cluster *the_cluster, struct star *stars, int index);
void calcAbsCoeffs (int filterSet);

#endif
