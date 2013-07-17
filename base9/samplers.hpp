#ifndef SAMPLERS_H
#define SAMPLERS_H

#include <random>

#include "densities.hpp"
#include "structures.hpp"

const double AGEPROPSTEPSIZE = 0.02;
const int AGE_WANDER         = 100;

// void propFieldStar (Star *inputStar);
// void propMass (Star *inputstar);
// void propMassRatio (Star *inputstar);
// void propClustParam (Cluster *clust, int TYPE);
// double gen_norm (double mean, double std_dev);
double sampleT (std::mt19937&, double var, double nu = DOF);
// double gamdev (double a);

#endif
