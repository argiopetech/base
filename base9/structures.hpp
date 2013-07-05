#ifndef STRUCT_H
#define STRUCT_H

#include <vector>

#include "constants.hpp"
#include "Cluster.hpp"
#include "Star.hpp"
#include "Chain.hpp"

struct globalIso
{
    globalIso()
        : globalIso(370, FILTS)
    {;}

    globalIso(int nEntries, int nFilts)
        : nEntries(nEntries),nFilts(nFilts)
    {
        assert(nEntries != 0 && nFilts != 0);

        mass.resize(nEntries);
        massNow.resize(nEntries);
        eep.resize(nEntries);

        mag.assign(nEntries, std::vector<double>(nFilts));
    }

    ~globalIso()
    {;}

    std::vector< std::vector<double> > mag;

    std::vector<int> eep;
    std::vector<double> mass;
    std::vector<double> massNow;

    int nEntries;
    int nFilts;

    double AgbTurnoffMass;
    double FeH;
    double logAge;
    double age;
    double Y;
    double z;
};

enum class blockType
{
    INITIAL_WANDER_BLOCK,
    MASS_RATIO_BLOCK,
    AGE_BLOCK,
    MODULUS_BLOCK,
    AGE_MODULUS_BLOCK,
    FEH_BLOCK,
    Y_BLOCK,
    TURN_FS_SAMPLING_OFF_BLOCK,
    FIELD_STAR_PRIOR_BLOCK,
    TURN_FS_SAMPLING_ON_BLOCK,
    FINAL_WANDER_BLOCK,
    END_BURN_IN_BLOCK,
    POST_BURN_IN_BLOCK
};

struct mcmcControl
{
    FILE *rData;
    FILE *wMass1File[2];
    FILE *wMass2File[2];
    FILE *wClusterFile[2];
    FILE *wDebugFile[2];
    FILE *rBetaInFile;
    FILE *wBetaOutFile;
    int fsSamplingOn;
    int sampleVarScale;
    int B0;                     // length of first "wandering" period
    int B2;                     // length of final "wandering" period
    int T;                      // length of decorrelation periods
    int increment;              // = ctrl.T / N_LSQ_VALUES
    int nIter;                  // number of post burn-in iterations
    int burnIter;                       // total number of burn-in iterations
    int runBurnIn;
    int outputBurnIn;
    int thin;
    double initialAge;
    double minMag;
    double maxMag;
    int iMag;
    int iStart;
    double filterPriorMin[FILTS];
    double filterPriorMax[FILTS];
    int useFilt[FILTS];
    enum blockType currentBlock;
};

struct block
{
    enum blockType type;
    void (*run) (int *iter, int nIter, Chain * mc, struct mcmcControl * ctrl);
    int nIter;
};

void swapGlobalEntries (struct globalIso &thisIso, int n, int useFilt[FILTS]);

#endif
