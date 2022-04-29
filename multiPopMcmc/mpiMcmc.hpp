#ifndef IFMRMCMC_H
#define IFMRMCMC_H

#include <array>
#include <string>
#include <fstream>

#include "Cluster.hpp"
#include "constants.hpp"
#include "Matrix.hpp"
#include "Star.hpp"

// WD evaluation constants
// For use in logPostStep
const double MIN_MASS1    = 0.15;

// The "global replacement structure" for mpiMcmc.
// This should eventually be deprecated, because it's awful.
struct ifmrMcmcControl
{
    int nIter;                  // number of post burn-in iterations
    int burnIter;               // total number of burn-in iterations
    int thin;
    std::array<double, NPARAMS> priorVar; // Rewrite the way this is handled... This is ugly.
};

#endif
