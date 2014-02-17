#ifndef SAMPLERS_H
#define SAMPLERS_H

#include <random>

#include "densities.hpp"
#include "structures.hpp"

const double AGEPROPSTEPSIZE = 0.02;
const int AGE_WANDER         = 100;

double sampleT (std::mt19937&, double var, double nu = DOF);

#endif
