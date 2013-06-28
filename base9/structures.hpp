#ifndef STRUCT_H
#define STRUCT_H

#include <vector>

#include "constants.hpp"
#include "Cluster.hpp"
#include "Star.hpp"

//Be careful adding new sample types.  Their order is important.
//There are a few places in the code that test for SAMPLE_TYPE > or <.
const int MASS              = -2;
const int AGE_DURING_WANDER = -1;  // age sampling during AGE_WANDER (defined in samplers.h)

const int AGE               =  0;  // age sampling
const int YYY               =  1;  // helium sampling
const int FEH               =  2;  // metallicity sampling
const int MOD               =  3;  // modulus sampling
const int ABS               =  4;  // absorption sampling;
const int IFMR_INTERCEPT    =  5;
const int IFMR_SLOPE        =  6;
const int IFMR_QUADCOEF     =  7;

struct globalIso
{
    double FeH;
    double logAge;
    double age;
    double Y;
    double z;
    double *mass;
    double *massNow;
    double **mag;                       //Dimensions are [NMASSES][NFILTS];
    double AgbTurnoffMass;
    int nEntries;
    int nFilts;
    int *eep;
};

struct chain
{
    Cluster clust;
    std::vector<Star> stars;
    double temperature;
    int acceptClust[NPARAMS];
    int rejectClust[NPARAMS];
    int *acceptMass;
    int *rejectMass;
    int *acceptMassRatio;
    int *rejectMassRatio;
    int *isFieldStarNow;
    int *isClusterStarNow;
};

enum blockType
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
    int verbose;
    int useFilt[FILTS];
    enum blockType currentBlock;
};

struct block
{
    enum blockType type;
    void (*run) (int *iter, int nIter, struct chain * mc, struct mcmcControl * ctrl);
    int nIter;
};


// Helper functions for star and cluster structures.
void readStar (FILE * pFile, Star *pStar);
void writeStar (FILE * pFile, Star *pStar);
double getMass1 (Star *pStar, Cluster *pCluster);
double getMass2 (Star *pStar, Cluster *pCluster);
void setMass1 (Star *pStar, Cluster *pCluster, double newMass);
void setMass2 (Star *pStar, Cluster *pCluster, double newMass);
void quickCopy (Star *pStarFrom, Star *pStarTo);
double getParameter (Cluster *pCluster, int TYPE);
void setParameter (Cluster *pCluster, int TYPE, double newValue);
void setFilterNames (int filterSet);
char *getFilterName (int index);
void allocateGlobalIso (struct globalIso *newIso);
void freeGlobalIso (struct globalIso *newIso);
void swapGlobalEntries (struct globalIso *thisIso, int n, int useFilt[FILTS]);

#endif
