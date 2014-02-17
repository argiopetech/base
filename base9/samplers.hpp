#ifndef SAMPLERS_H
#define SAMPLERS_H

#include <random>

#include "constants.hpp"

double sampleT (std::mt19937&, double var, double nu = DOF);

#endif
