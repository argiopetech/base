#ifndef IFMRMCMC_H
#define IFMRMCMC_H

#include <array>
#include <string>
#include <fstream>

#include "Cluster.hpp"
#include "constants.hpp"
#include "Matrix.hpp"
#include "Star.hpp"

const double MIN_MASS1    = 0.15;
const int N_WD_MASS1      = 8000;  /* evaluate white dwarfs on a finer grid? */

struct ifmrMcmcControl
{
    std::string clusterFilename;

    int nIter;                  // number of post burn-in iterations
    int burnIter;               // total number of burn-in iterations
    int thin;
    std::array<double, NPARAMS> priorVar; // Rewrite the way this is handled... This is ugly.
    Matrix<double, NPARAMS, NPARAMS> propMatrix;
};

#endif
