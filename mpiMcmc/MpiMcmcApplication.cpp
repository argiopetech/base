#include <iostream>

#include <boost/format.hpp>

#include "marg.hpp"
#include "mpiMcmc.hpp"
#include "mpiFuncs.hpp"
#include "MpiMcmcApplication.hpp"
#include "Settings.hpp"
#include "samplers.cpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::cout;
using std::cerr;
using std::endl;

int MpiMcmcApplication::run()
{
    int increment;

    double logPostCurr;
    double logPostProp;
    double fsLike;

    Chain mc;
    struct ifmrMcmcControl ctrl;
    Cluster propClust;

    array<double, 2> ltau;

    array<double, FILTS> filterPriorMin;
    array<double, FILTS> filterPriorMax;

    Matrix<double, NPARAMS, nSave> params;

    std::vector<int> filters;

    increment = settings.mpiMcmc.burnIter / (2 * nSave);

    initIfmrMcmcControl (mc.clust, ctrl, evoModels, settings);

    /* Initialize filter prior mins and maxes */
    filterPriorMin.fill(1000);
    filterPriorMax.fill(-1000);

    readCmdData (mc.stars, ctrl, evoModels, filters, filterPriorMin, filterPriorMax, settings);

    initChain (mc, evoModels, ltau, filters);

    if (mc.stars.size() > 1)
    {
        double logFieldStarLikelihood = 0.0;

        for (decltype(filters.size()) filt = 0; filt < filters.size(); filt++)
        {
            logFieldStarLikelihood -= log (filterPriorMax[filt] - filterPriorMin[filt]);
        }
        fsLike = exp (logFieldStarLikelihood);
    }
    else
    {
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
        if (iteration < ctrl.burnIter / 2)
        {
            propClust = propClustBigSteps (mc.clust, ctrl);
        }
        else
        {
            propClust = propClustIndep (mc.clust, ctrl);
        }

        try
        {
            logPostProp = logPostStep (mc, propClust, fsLike, filters, filterPriorMin, filterPriorMax);
        }
        catch(InvalidCluster &e)
        {
            logPostProp = -HUGE_VAL;
        }

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp))
        {
            mc.clust = propClust;
            logPostCurr = logPostProp;
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
                        params.at(p).at((iteration - ctrl.burnIter / 2) / increment) = mc.clust.getParam(p);
                    }
                }
            }
        }

        /* Write output */
        for (int p = 0; p < NPARAMS; p++)
        {
            if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)
            {
                ctrl.burninFile << boost::format("%10.6f ") % mc.clust.getParam(p);
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
        propClust = propClustCorrelated (mc.clust, ctrl);

        try
        {
            logPostProp = logPostStep (mc, propClust, fsLike, filters, filterPriorMin, filterPriorMax);
        }
        catch(InvalidCluster &e)
        {
            logPostProp = -HUGE_VAL;
        }

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp))
        {
            mc.clust = propClust;
            logPostCurr = logPostProp;
        }

        if (iteration % ctrl.thin == 0)
        {
            for (int p = 0; p < NPARAMS; p++)
            {
                if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)
                {
                    ctrl.resFile << boost::format("%10.6f ") % mc.clust.getParam(p);
                }
            }
            ctrl.resFile << boost::format("%10.6f") % logPostCurr << endl;
        }
    }

    ctrl.resFile.close();

    cout << "\nAcceptance ratio: " << acceptanceRatio() << endl;

    return 0;
}


Cluster MpiMcmcApplication::propClustBigSteps (Cluster clust, struct ifmrMcmcControl const &ctrl)
{
    return propClustIndep(clust, ctrl, 25.0);
}

Cluster MpiMcmcApplication::propClustIndep (Cluster clust, struct ifmrMcmcControl const &ctrl, double scale)
{
    /* DOF defined in densities.h */
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            clust.setParam(p, clust.getParam(p) + sampleT (gen, scale * clust.stepSize.at(p) * clust.stepSize.at(p)));
        }
    }

    return clust;
}

Cluster MpiMcmcApplication::propClustCorrelated (Cluster clust, struct ifmrMcmcControl const &ctrl)
{
    /* DOF defined in densities.h */
    array<double, NPARAMS> indepProps;
    array<double, NPARAMS> corrProps;

    indepProps.fill(0.0);
    corrProps.fill(0.0);

    int p, k;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            indepProps.at(p) = sampleT (gen, 1.0);
        }
    }
    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            for (k = 0; k <= p; k++)
            {                           /* propMatrix is lower diagonal */
                if (ctrl.priorVar.at(k) > EPSILON)
                {
                    corrProps.at(p) += ctrl.propMatrix.at(p).at(k) * indepProps.at(k);
                }
            }
            clust.setParam(p, clust.getParam(p) + corrProps.at(p));
        }
    }

    return clust;
}

double MpiMcmcApplication::logPostStep(Chain &mc, Cluster &propClust, double fsLike, const vector<int> &filters, std::array<double, FILTS> &filterPriorMin, std::array<double, FILTS> &filterPriorMax)
{
    mutex logPostMutex;
    double logPostProp;

    logPostProp = logPriorClust (propClust, evoModels);

    auto stars = mc.stars;

    propClust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, propClust.feh, propClust.yyy, propClust.age);    // determine AGBt ZAMS mass, to find evol state

    /* loop over assigned stars */
    pool.parallelFor(mc.stars.size(), [=,&logPostMutex,&logPostProp](int i)
    {
        double postClusterStar = 0.0;

        array<double, FILTS> globalMags;
        array<double, 2> ltau;

        /* loop over all (mass1, mass ratio) pairs */
        if (stars.at(i).status[0] == WD)
        {
            postClusterStar = 0.0;

            for (int j = 0; j < N_WD_MASS1; j++)
            {
                double tmpLogPost;
                Star wd(stars.at(i));
                wd.isFieldStar = 0;
                wd.U = wdGridMass(j);
                wd.massRatio = 0.0;

                evolve (propClust, evoModels, globalMags, filters, wd, ltau);

                try
                {
                    tmpLogPost = logPost1Star (wd, propClust, evoModels, filterPriorMin, filterPriorMax);
                    tmpLogPost += log ((mc.clust.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1);

                    postClusterStar = postClusterStar +  exp (tmpLogPost);
                }
                catch ( WDBoundsError &e )
                {
                    cerr << e.what() << endl;
                }
            }
        }
        else
        {
            /* marginalize over isochrone */
            postClusterStar = margEvolveWithBinary (propClust, stars.at(i), evoModels, filters, ltau, globalMags, filterPriorMin, filterPriorMax);
        }
            
        postClusterStar = postClusterStar * stars.at(i).clustStarPriorDens;


        /* marginalize over field star status */
        std::lock_guard<mutex> lk(logPostMutex);
        logPostProp += log ((1.0 - stars.at(i).clustStarPriorDens) * fsLike + postClusterStar);
    });

    return logPostProp;
}

double MpiMcmcApplication::wdGridMass (int point) const
{
    // I think this only gets calculated once, but I can't quite be sure...
    // It may not even matter, but every little bit helps, eh? And hey, preoptimization...
    static const double dWdMass1 = (settings.whiteDwarf.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1;

    return MIN_MASS1 + (dWdMass1 * point);
}
