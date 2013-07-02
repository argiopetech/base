#include <gtest/gtest.h>

#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/format.hpp>

#include "evolve.hpp"
#include "mt19937ar.hpp"
#include "Settings.hpp"
#include "loadModels.hpp"
#include "FilterSet.hpp"
#include "mpi_things.hpp"

#include "../mpiFuncs.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::array;
using std::string;
using std::ofstream;

extern int numFilts;

Model yamlChunk();
double run1step(void);

TEST(mpiMcmc, oneStep)
{
    EXPECT_NEAR(-448.007816, run1step(), 0.0001);
}

void initMassGrids (array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &msMass1Grid, array<double, N_MS_MASS1 * N_MS_MASS_RATIO> &msMassRatioGrid, array<double, N_WD_MASS1> &wdMass1Grid, Chain const &mc);
void initChain (Chain &mc, const struct ifmrMcmcControl &ctrl, Model &evoModels, array<double, 2>&);
void readCmdData (Chain &mc, struct ifmrMcmcControl &ctrl, Model &evoModels);
void initIfmrMcmcControl (Chain &mc, struct ifmrMcmcControl &ctrl, Model &evoModels);
void initStepSizes (Cluster &clust);
void printHeader (ofstream &file, array<double, NPARAMS> const &priors);
void propClustCorrelated (Cluster &clust, struct ifmrMcmcControl const &ctrl);
void propClustIndep (Cluster &clust, struct ifmrMcmcControl const &ctrl);
void propClustBigSteps (Cluster &clust, struct ifmrMcmcControl const &ctrl);

Model yamlChunk()
{
    settings.fromYaml ("/home/elliot/Projects/stellar_evolution/test/hyades2/base9.yaml");

    return makeModel(settings);
}


double run1step()
{
    double logPostProp;
    double fsLike;

    Chain mc;
    struct ifmrMcmcControl ctrl;
    Cluster propClust;

    array<double, 2> ltau;
    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMass1Grid;
    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMassRatioGrid;
    array<double, N_WD_MASS1> wdMass1Grid;

    settings.fromYaml ("/home/elliot/Projects/stellar_evolution/test/hyades2/base9.yaml");

    Model evoModels = makeModel(settings);

    settings.files.output = "/home/elliot/Projects/stellar_evolution/test/hyades2/hyades/Mcmc";
    settings.files.phot   = "/home/elliot/Projects/stellar_evolution/test/hyades2/hyades/Hyades.UBV.phot";
    settings.files.models = "/home/elliot/Projects/stellar_evolution/test/hyades2/models/";

    initStepSizes (mc.clust);

    initIfmrMcmcControl (mc, ctrl, evoModels);

    for (int p = 0; p < NPARAMS; p++)
    {
        mc.clust.priorVar[p] = ctrl.priorVar[p];
        mc.clust.priorMean[p] = ctrl.priorMean[p];
    }

    readCmdData (mc, ctrl, evoModels);

    evoModels.numFilts = ctrl.numFilts;
    numFilts = ctrl.numFilts;

    initChain (mc, ctrl, evoModels, ltau);

    for (int i = 0; i < mc.clust.nStars; i++)
    {
        mc.stars[i].isFieldStar = 0;
        mc.stars[i].boundsFlag = 0;
    }

    initMassGrids (msMass1Grid, msMassRatioGrid, wdMass1Grid, mc);

    double logFieldStarLikelihood = 0.0;

    for (int filt = 0; filt < ctrl.numFilts; filt++)
    {
        logFieldStarLikelihood -= log (ctrl.filterPriorMax[filt] - ctrl.filterPriorMin[filt]);
    }
    fsLike = exp (logFieldStarLikelihood);

    printHeader (ctrl.burninFile, ctrl.priorVar);

    propClust = mc.clust;

    propClustBigSteps (propClust, ctrl);

    if (ctrl.priorVar[ABS] > EPSILON)
    {
        propClust.parameter[ABS] = fabs (propClust.parameter[ABS]);
    }
    if (ctrl.priorVar[IFMR_SLOPE] > EPSILON)
    {
        propClust.parameter[IFMR_SLOPE] = fabs (propClust.parameter[IFMR_SLOPE]);
    }

    logPostProp = logPostStep (mc, evoModels, wdMass1Grid, propClust, fsLike, ltau);

    return logPostProp;
}
