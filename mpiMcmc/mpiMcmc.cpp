#include <boost/format.hpp>

#include "mpiMcmc.hpp"
#include "mpiFuncs.hpp"
#include "Settings.hpp"

using std::array;
using std::cout;
using std::cerr;
using std::endl;

int main (int argc, char *argv[])
{
    int accept = 0, reject = 0;
    int increment;

    double logPostCurr;
    double logPostProp;
    double fsLike;

    Chain mc;
    struct ifmrMcmcControl ctrl;
    Cluster propClust;

    array<double, 2> ltau;
    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMass1Grid;
    array<double, N_MS_MASS1 * N_MS_MASS_RATIO> msMassRatioGrid;
    array<double, N_WD_MASS1> wdMass1Grid;

    Matrix<double, NPARAMS, nSave> params;

    Settings settings;

    // Setup settings
    {
        settings.fromCLI (argc, argv);
        if (!settings.files.config.empty())
        {
            settings.fromYaml (settings.files.config);
        }
        else
        {
            settings.fromYaml ("base9.yaml");
        }

        settings.fromCLI (argc, argv);
    }

    const Model evoModels = makeModel(settings);

    increment = settings.mpiMcmc.burnIter / (2 * nSave);

    initIfmrMcmcControl (mc, ctrl, evoModels, settings);

    for (int p = 0; p < NPARAMS; p++)
    {
        mc.clust.priorVar[p] = ctrl.priorVar[p];
        mc.clust.priorMean[p] = ctrl.priorMean[p];
    }

    readCmdData (mc, ctrl, evoModels);

    initChain (mc, ctrl, evoModels, ltau);

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

        logPostProp = logPostStep (mc, evoModels, wdMass1Grid, propClust, fsLike, ltau);

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp, ltau))
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
    }

    ctrl.burninFile.close();

    make_cholesky_decomp(ctrl, params);

    // Main run
    ctrl.resFile.open(ctrl.clusterFilename);
    if (!ctrl.resFile)
    {
        cerr << "***Error: File " << ctrl.clusterFilename << " was not available for writing.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    printHeader (ctrl.resFile, ctrl.priorVar);

    for (int iteration = 0; iteration < ctrl.nIter * ctrl.thin; iteration++)
    {
        propClust = mc.clust;
        propClustCorrelated (propClust, ctrl);

        logPostProp = logPostStep (mc, evoModels, wdMass1Grid, propClust, fsLike, ltau);

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp, ltau))
        {
            mc.clust = propClust;
            logPostCurr = logPostProp;
            accept++;
        }
        else
        {
            reject++;
        }

        if (iteration % ctrl.thin == 0)
        {
            for (int p = 0; p < NPARAMS; p++)
            {
                if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)
                {
                    ctrl.resFile << boost::format("%10.6f ") % mc.clust.parameter[p];
                }
            }
            ctrl.resFile << boost::format("%10.6f") % logPostCurr << endl;
        }
    }

    ctrl.resFile.close();

    cout << "\nAcceptance ratio: " << (double) accept / (accept + reject) << endl;

    return 0;
}
