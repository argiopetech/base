#ifndef SAMPLERS_H
#define SAMPLERS_H

#include "structures.hpp"

const double AGEPROPSTEPSIZE = 0.02;
const int AGE_WANDER         = 100;

void propFieldStar (struct star *inputStar);
void propMass (struct star *inputstar);
void propMassRatio (struct star *inputstar);
void propClustParam (struct cluster *clust, int TYPE);
double gen_norm (double mean, double std_dev);
double sampleT (double var, double nu);
double gamdev (double a);

#endif
