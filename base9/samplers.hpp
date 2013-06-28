#ifndef SAMPLERS_H
#define SAMPLERS_H

#include "structures.hpp"

const double AGEPROPSTEPSIZE = 0.02;
const int AGE_WANDER         = 100;

void propFieldStar (Star *inputStar);
void propMass (Star *inputstar);
void propMassRatio (Star *inputstar);
void propClustParam (Cluster *clust, int TYPE);
double gen_norm (double mean, double std_dev);
double sampleT (double var, double nu);
double gamdev (double a);

#endif
