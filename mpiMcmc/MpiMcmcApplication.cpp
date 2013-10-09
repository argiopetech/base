#include <iostream>

#include <boost/format.hpp>

#include "Cluster.hpp"
#include "marg.hpp"
#include "mpiMcmc.hpp"
#include "mpiFuncs.hpp"
#include "MpiMcmcApplication.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "samplers.cpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::cout;
using std::cerr;
using std::endl;

MpiMcmcApplication::MpiMcmcApplication(Settings &s)
        : McmcApplication(s.seed), evoModels(makeModel(s)), settings(s), pool(s.threads)
{
    ctrl.priorVar.fill(0);

    clust.feh = clust.priorMean[FEH] = settings.cluster.Fe_H;
    ctrl.priorVar[FEH] = settings.cluster.sigma.Fe_H;

    clust.mod = clust.priorMean[MOD] = settings.cluster.distMod;
    ctrl.priorVar[MOD] = settings.cluster.sigma.distMod;

    clust.abs = clust.priorMean[ABS] = fabs(s.cluster.Av);
    ctrl.priorVar[ABS] = settings.cluster.sigma.Av;

    clust.age = clust.priorMean[AGE] = settings.cluster.logClusAge;
    ctrl.priorVar[AGE] = 1.0;

    if (s.whiteDwarf.wdModel == WdModel::MONTGOMERY)
    {
        clust.carbonicity = clust.priorMean[CARBONICITY] = settings.cluster.carbonicity;
        ctrl.priorVar[CARBONICITY] = settings.cluster.sigma.carbonicity;
    }
    else
    {
        clust.carbonicity = clust.priorMean[CARBONICITY] = 0.0;
        ctrl.priorVar[CARBONICITY] = 0.0;
    }

    if (s.mainSequence.msRgbModel == MsModel::CHABHELIUM)
    {
        clust.yyy = clust.priorMean[YYY] = settings.cluster.Y;
        ctrl.priorVar[YYY] = settings.cluster.sigma.Y;
    }
    else
    {
        clust.yyy = clust.priorMean[YYY] = 0.0;
        ctrl.priorVar[YYY] = 0.0;
    }


    if (evoModels.IFMR <= 3)
    {
        ctrl.priorVar[IFMR_SLOPE] = 0.0;
        ctrl.priorVar[IFMR_INTERCEPT] = 0.0;
        ctrl.priorVar[IFMR_QUADCOEF] = 0.0;
    }
    else if (evoModels.IFMR <= 8)
    {
        ctrl.priorVar[IFMR_SLOPE] = 1.0;
        ctrl.priorVar[IFMR_INTERCEPT] = 1.0;
        ctrl.priorVar[IFMR_QUADCOEF] = 0.0;
    }
    else
    {
        ctrl.priorVar[IFMR_SLOPE] = 1.0;
        ctrl.priorVar[IFMR_INTERCEPT] = 1.0;
        ctrl.priorVar[IFMR_QUADCOEF] = 1.0;
    }

    /* set starting values for IFMR parameters */
    clust.ifmrSlope = clust.priorMean[IFMR_SLOPE] = 0.08;
    clust.ifmrIntercept = clust.priorMean[IFMR_INTERCEPT] = 0.65;

    if (evoModels.IFMR <= 10)
        clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.0001;
    else
        clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.08;

    for (auto &var : ctrl.priorVar)
    {
        if (var < 0.0)
        {
            var = 0.0;
        }
        else
        {
            var = var * var;
        }
    }

    /* read burnIter and nIter */
    {
        ctrl.burnIter = settings.mpiMcmc.burnIter;
        int nParamsUsed = 0;

        for (int p = 0; p < NPARAMS; p++)
        {
            if (ctrl.priorVar.at(p) > EPSILON)
            {
                nParamsUsed++;
            }
        }

        if (ctrl.burnIter < 2 * (nParamsUsed + 1))
        {
            ctrl.burnIter = 2 * (nParamsUsed + 1);
            cerr << "burnIter below minimum allowable size. Increasing to " << ctrl.burnIter << endl;
        }

        if ((ctrl.burnIter / 2) < 100)
        {
            nSave = ctrl.burnIter / 2;
        }
    }

    ctrl.nIter = settings.mpiMcmc.maxIter;
    ctrl.thin = settings.mpiMcmc.thin;

    /* open files for reading (data) and writing */
    string filename;

    filename = settings.files.phot;
    ctrl.rData.open(filename);
    if (!ctrl.rData)
    {
        cerr << "***Error: Photometry file " << filename << " was not found.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    if (s.cluster.index < 0 || settings.cluster.index > FILTS)
    {
        cerr << "***Error: " << settings.cluster.index << " not a valid magnitude index.  Choose 0, 1,or 2.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    ctrl.clusterFilename = settings.files.output + ".res";

    std::copy(ctrl.priorVar.begin(), ctrl.priorVar.end(), clust.priorVar.begin());
}


int MpiMcmcApplication::run()
{
   int increment;

    double logPostCurr;
    double logPostProp;
    double fsLike;

    Cluster propClust;

    array<double, 2> ltau;

    array<double, FILTS> filterPriorMin;
    array<double, FILTS> filterPriorMax;

    MVatrix<double, NPARAMS> params; // Must be initialized after nSave has been set.

    std::vector<int> filters;

    params.fill(vector<double>(nSave, 0.0));

    increment = settings.mpiMcmc.burnIter / (2 * nSave);

    /* Initialize filter prior mins and maxes */
    filterPriorMin.fill(1000);
    filterPriorMax.fill(-1000);

    stars = readCmdData (ctrl, evoModels, filters, filterPriorMin, filterPriorMax, settings);

    // Begin initChain
    {
        array<double, FILTS> globalMags;

        for (auto star : stars)
        {
            star.clustStarProposalDens = star.clustStarPriorDens;   // Use prior prob of being clus star
            star.UStepSize = 0.001; // within factor of ~2 for most main sequence stars
            star.massRatioStepSize = 0.001;

            // find photometry for initial values of currentClust and mc.stars
            clust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, clust.feh, clust.yyy, clust.age);    // determine AGBt ZAMS mass, to find evol state
            evolve (clust, evoModels, globalMags, filters, star, ltau);

            if (star.status[0] == WD)
            {
                star.UStepSize = 0.05;      // use larger initial step size for white dwarfs
                star.massRatio = 0.0;
            }
        }
    }
    // end initChain

    if (stars.size() > 1)
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

    /* set current log posterior to minimum double value */
    /* will cause random starting value */
    logPostCurr = -HUGE_VAL; //std::numeric_limits<double>::min();

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
        if (settings.mpiMcmc.bigStepBurnin || (iteration < ctrl.burnIter / 2))
        {
            propClust = propClustBigSteps (clust, ctrl);
        }
        else
        {
            propClust = propClustIndep (clust, ctrl);
        }

        try
        {
            logPostProp = logPostStep (stars, propClust, fsLike, filters, filterPriorMin, filterPriorMax);
        }
        catch(InvalidCluster &e)
        {
            logPostProp = -HUGE_VAL;
        }

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp))
        {
            clust = propClust;
            logPostCurr = logPostProp;
        }

        /* save draws to estimate covariance matrix for more efficient Metropolis */
        if ((iteration >= (ctrl.burnIter - (nSave * increment))) && (iteration < ctrl.burnIter))
        {
            if (iteration % increment == 0)
            {
                /* save draws */
                for (int p = 0; p < NPARAMS; p++)
                {
                    if (ctrl.priorVar[p] > EPSILON)
                    {
                        params.at(p).at((iteration - (ctrl.burnIter - (nSave * increment))) / increment) = clust.getParam(p);
                    }
                }
            }
        }

        /* Write output */
        for (int p = 0; p < NPARAMS; p++)
        {
            if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)
            {
                ctrl.burninFile << boost::format("%10.6f ") % clust.getParam(p);
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
        propClust = propClustCorrelated (clust, ctrl);

        try
        {
            logPostProp = logPostStep (stars, propClust, fsLike, filters, filterPriorMin, filterPriorMax);
        }
        catch(InvalidCluster &e)
        {
            logPostProp = -HUGE_VAL;
        }

        /* accept/reject */
        if (acceptClustMarg (logPostCurr, logPostProp))
        {
            clust = propClust;
            logPostCurr = logPostProp;
        }

        if (iteration % ctrl.thin == 0)
        {
            for (int p = 0; p < NPARAMS; p++)
            {
                if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)
                {
                    ctrl.resFile << boost::format("%10.6f ") % clust.getParam(p);
                }
            }
            ctrl.resFile << boost::format("%10.6f") % logPostCurr << endl;
        }
    }

    ctrl.resFile.close();

    cout << "\nAcceptance ratio: " << acceptanceRatio() << endl;

    return 0;
}


Cluster MpiMcmcApplication::propClustBigSteps (const Cluster &clust, struct ifmrMcmcControl const &ctrl)
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

double scaleFactor (double acceptanceRatio)
{
    double factor = 1.0;

    if (acceptanceRatio > 0.9)
    {
        factor = 2.0;
    }
    else if (acceptanceRatio > 0.7)
    {
        factor = 1.8;
    }
    else if (acceptanceRatio > 0.5)
    {
        factor = 1.5;
    }
    else if (acceptanceRatio > 0.4)
    {
        factor = 1.2;
    }
    else if (acceptanceRatio < 0.2)
    {
        factor = 1 / 1.5;
    }
    else if (acceptanceRatio < 0.15)
    {
        factor = 1 / 1.8;
    }
    else if (acceptanceRatio < 0.05)
    {
        factor = 0.5;
    }

    return factor;
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

double MpiMcmcApplication::logPostStep(const vector<Star> &stars, Cluster &propClust, double fsLike, const vector<int> &filters, std::array<double, FILTS> &filterPriorMin, std::array<double, FILTS> &filterPriorMax)
{
    mutex logPostMutex;
    double logPostProp;

    logPostProp = logPriorClust (propClust, evoModels);

    propClust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, propClust.feh, propClust.yyy, propClust.age);    // determine AGBt ZAMS mass, to find evol state

    /* loop over assigned stars */
    pool.parallelFor(stars.size(), [=,&logPostMutex,&logPostProp](int i)
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
                    tmpLogPost += log ((propClust.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1);

                    postClusterStar +=  exp (tmpLogPost);
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

        postClusterStar *= stars.at(i).clustStarPriorDens;

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
