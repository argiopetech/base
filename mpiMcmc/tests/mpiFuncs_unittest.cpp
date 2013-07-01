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
void initChain (Chain &mc, const struct ifmrMcmcControl &ctrl, Model &evoModels);
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
    int accept = 0, reject = 0;
    int increment;

    double logPostCurr;
    double logPostProp;
    double fsLike;

    Chain mc;
    struct ifmrMcmcControl ctrl;
    Cluster propClust;

    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMass1Grid;
    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMassRatioGrid;
    array<double, N_WD_MASS1> wdMass1Grid;

    Matrix<double, NPARAMS, nSave> params;

    settings.fromYaml ("/home/elliot/Projects/stellar_evolution/test/hyades2/base9.yaml");

    Model evoModels = makeModel(settings);

    settings.files.output = "/home/elliot/Projects/stellar_evolution/test/hyades2/hyades/Mcmc";
    settings.files.phot = "/home/elliot/Projects/stellar_evolution/test/hyades2/hyades/Hyades.UBV.phot";
    settings.files.models = "/home/elliot/Projects/stellar_evolution/test/hyades2/models/";

    increment = settings.mpiMcmc.burnIter / (2 * nSave);

    initStepSizes (mc.clust);

    initIfmrMcmcControl (mc, ctrl, evoModels);

    evoModels.WDatm = BERGERON;

    for (int p = 0; p < NPARAMS; p++)
    {
        mc.clust.priorVar[p] = ctrl.priorVar[p];
        mc.clust.priorMean[p] = ctrl.priorMean[p];
    }

    readCmdData (mc, ctrl, evoModels);

    evoModels.numFilts = ctrl.numFilts;
    numFilts = ctrl.numFilts;

    initChain (mc, ctrl, evoModels);

    for (int i = 0; i < mc.clust.nStars; i++)
    {
        mc.stars[i].isFieldStar = 0;
        mc.stars[i].boundsFlag = 0;
    }

    initMassGrids (msMass1Grid, msMassRatioGrid, wdMass1Grid, mc);

    double logFieldStarLikelihood = 0.0;

    if (mc.clust.nStars > 1)
    {
        for (int filt = 0; filt < ctrl.numFilts; filt++)
        {
            logFieldStarLikelihood -= log (ctrl.filterPriorMax[filt] - ctrl.filterPriorMin[filt]);
        }
        fsLike = exp (logFieldStarLikelihood);
    }
    else
    {
        logFieldStarLikelihood = -HUGE_VAL;
        fsLike = 0;
    }

    cout << "Bayesian analysis of stellar evolution" << endl;

    /* set current log posterior to -HUGE_VAL */
    /* will cause random starting value */
    logPostCurr = -HUGE_VAL;

    // Run Burnin
    ctrl.burninFile.open(ctrl.clusterFilename + ".burnin");
    if (!ctrl.burninFile)
    {
        cerr << "***Error: File " << ctrl.clusterFilename << " was not available for writing.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    printHeader (ctrl.burninFile, ctrl.priorVar);

    for (int iteration = 0; iteration < ctrl.burnIter; iteration++)
    {
        propClust = mc.clust;

        if (iteration < ctrl.burnIter / 2)
        {
            propClustBigSteps (propClust, ctrl);
        }
        else
        {
            propClustIndep (propClust, ctrl);
        }

        if (ctrl.priorVar[ABS] > EPSILON)
        {
            propClust.parameter[ABS] = fabs (propClust.parameter[ABS]);
        }
        if (ctrl.priorVar[IFMR_SLOPE] > EPSILON)
        {
            propClust.parameter[IFMR_SLOPE] = fabs (propClust.parameter[IFMR_SLOPE]);
        }

        logPostProp = logPostStep (mc, evoModels, wdMass1Grid, propClust, fsLike);

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp))
        {
            mc.clust = propClust;
            logPostCurr = logPostProp;
            accept++;
        }
        else
        {
            reject++;
        }
        /* save draws to estimate covariance matrix for more efficient Metropolis */
        if (iteration >= ctrl.burnIter / 2 && iteration < ctrl.burnIter)
        {
            if (iteration % increment == 0)
            {
                /* save draws */
                for (int p = 0; p < NPARAMS; p++)
                {
                    if (ctrl.priorVar[p] > EPSILON)
                    {
                        params.at(p).at((iteration - ctrl.burnIter / 2) / increment) = mc.clust.parameter[p];
                    }
                }
            }
        }

        /* Write output */
        for (int p = 0; p < NPARAMS; p++)
        {
            if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)
            {
                ctrl.burninFile << boost::format("%10.6f ") % mc.clust.parameter[p];
            }
        }

        ctrl.burninFile << boost::format("%10.6f") % logPostCurr << endl;
        return logPostCurr;
    }
}
