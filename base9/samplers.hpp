/* samplers.h */

#ifdef SAMPLERS_H
/* the file has been included already */
#else

#define SAMPLERS_H

#define AGEPROPSTEPSIZE         0.02
#define AGE_WANDER               100

#include "structures.hpp"

void propFieldStar (struct star *inputStar);
void propMass (struct star *inputstar);
void propMassRatio (struct star *inputstar);
void propClustParam (struct cluster *clust, int TYPE);
double gen_norm (double mean, double std_dev);
double sampleT (double var, double nu);
double gamdev (double a);

#endif
