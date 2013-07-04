#ifndef IFMRMCMC_H
#define IFMRMCMC_H

#include <array>
#include <string>
#include <fstream>

#include "evolve.hpp"

const double MIN_MASS1    = 0.15;
const int N_MS_MASS1      = 500;
const int N_MS_MASS_RATIO = 20;
const int N_WD_MASS1      = 8000;  /* evaluate white dwarfs on a finer grid? */
const int CLUS_FILE = 10;
const int ALLOC_CHUNK = 5;
const int nSave = 10;             /*changed from 100 to 10 */

struct obsStar
{
    double obsPhot[FILTS];
    double variance[FILTS];
    double clustStarPriorDens;  /* cluster membership prior probability */
};

struct ifmrMcmcControl
{
    std::ifstream rData;
    std::ofstream resFile;
    std::ofstream burninFile;
    std::string clusterFilename;

    int fsSamplingOn;
    int sampleVarScale;
    int nIter;                  // number of post burn-in iterations
    int burnIter;               // total number of burn-in iterations
    int thin;
    int modelSet;
    std::array<double, NPARAMS> priorVar;
    double initialAge;
    double minMag;
    double maxMag;
    int iMag;
    int iStart;
    double filterPriorMin[FILTS];
    double filterPriorMax[FILTS];
    int useFilt[FILTS];
    int numFilts;
    double propMatrix[NPARAMS][NPARAMS];
};

#endif
