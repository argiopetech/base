#if defined( EVOLVE_H )
// the file has been included already
#else
#define EVOLVE_H

#include "structures.hpp"

#define VERSION            4.8  // last update: 25aug09

void evolve (struct cluster *the_cluster, struct star *stars, int index);
void calcAbsCoeffs (int filterSet);

#endif
