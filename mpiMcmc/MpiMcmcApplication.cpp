#include <functional>
#include <iostream>

#include <boost/format.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#include "Chain.hpp"
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
using std::flush;
using std::function;
using std::mutex;
using std::vector;

using namespace std::placeholders;

MpiMcmcApplication::MpiMcmcApplication(Settings &s)
    : evoModels(makeModel(s)), settings(s), gen(uint32_t(s.seed * uint32_t(2654435761))), pool(s.threads) // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)
{
    ctrl.priorVar.fill(0);

    mainClust.feh = settings.cluster.starting.Fe_H;
    mainClust.priorMean[FEH] = clust.feh = clust.priorMean[FEH] = settings.cluster.Fe_H;
    ctrl.priorVar[FEH] = settings.cluster.sigma.Fe_H;

    mainClust.mod = settings.cluster.starting.distMod;
    mainClust.priorMean[MOD] = clust.mod = clust.priorMean[MOD] = settings.cluster.distMod;
    ctrl.priorVar[MOD] = settings.cluster.sigma.distMod;

    mainClust.abs = settings.cluster.starting.Av;
    mainClust.priorMean[ABS] = clust.abs = clust.priorMean[ABS] = fabs(settings.cluster.Av);
    ctrl.priorVar[ABS] = settings.cluster.sigma.Av;

    mainClust.age = mainClust.priorMean[AGE] = clust.age = clust.priorMean[AGE] = settings.cluster.logClusAge;
    ctrl.priorVar[AGE] = 1.0;

    if (settings.whiteDwarf.wdModel == WdModel::MONTGOMERY)
    {
        mainClust.carbonicity = settings.cluster.starting.carbonicity;
        mainClust.priorMean[CARBONICITY] = clust.carbonicity = clust.priorMean[CARBONICITY] = settings.cluster.carbonicity;
        ctrl.priorVar[CARBONICITY] = settings.cluster.sigma.carbonicity;
    }
    else
    {
        mainClust.carbonicity = mainClust.priorMean[CARBONICITY] = clust.carbonicity = clust.priorMean[CARBONICITY] = 0.0;
        ctrl.priorVar[CARBONICITY] = 0.0;
    }

    if (settings.mainSequence.msRgbModel == MsModel::CHABHELIUM)
    {
        mainClust.yyy = settings.cluster.starting.Y;
        mainClust.priorMean[YYY] = clust.yyy = clust.priorMean[YYY] = settings.cluster.Y;
        ctrl.priorVar[YYY] = settings.cluster.sigma.Y;
    }
    else
    {
        mainClust.yyy = mainClust.priorMean[YYY] = clust.yyy = clust.priorMean[YYY] = 0.0;
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
    mainClust.ifmrSlope = mainClust.priorMean[IFMR_SLOPE] =clust.ifmrSlope = clust.priorMean[IFMR_SLOPE] = 0.08;
    mainClust.ifmrIntercept = mainClust.priorMean[IFMR_INTERCEPT] = clust.ifmrIntercept = clust.priorMean[IFMR_INTERCEPT] = 0.65;

    if (evoModels.IFMR <= 10)
        mainClust.ifmrQuadCoef = mainClust.priorMean[IFMR_QUADCOEF] = clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.0001;
    else
        mainClust.ifmrQuadCoef = mainClust.priorMean[IFMR_QUADCOEF] = clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.08;

    clust.setM_wd_up(settings.whiteDwarf.M_wd_up);
    mainClust.setM_wd_up(settings.whiteDwarf.M_wd_up);

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
    }

    ctrl.nIter = settings.mpiMcmc.maxIter;
    ctrl.thin = settings.mpiMcmc.thin;

    ctrl.clusterFilename = settings.files.output + ".res";

    std::copy(ctrl.priorVar.begin(), ctrl.priorVar.end(), clust.priorVar.begin());
    std::copy(ctrl.priorVar.begin(), ctrl.priorVar.end(), mainClust.priorVar.begin());
}


int MpiMcmcApplication::run()
{
    double fsLike;

    array<double, NPARAMS> stepSize;
    std::copy(settings.mpiMcmc.stepSize.begin(), settings.mpiMcmc.stepSize.end(), stepSize.begin());

    {
        // Read photometry and calculate fsLike
        vector<double> filterPriorMin;
        vector<double> filterPriorMax;

        // open files for reading (data) and writing
        // rData implcitly relies on going out of scope to close the photometry file
        // This is awful, but pretty (since this code is, at time of writing, in restricted, anonymous scope
        std::ifstream rData(settings.files.phot);

        if (!rData)
        {
            cerr << "***Error: Photometry file " << settings.files.phot << " was not found.***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (-1);
        }

        auto ret = base::utility::readPhotometry (rData, filterPriorMin, filterPriorMax, settings);
        auto filterNames = ret.first;
        systems = ret.second;

        evoModels.restrictFilters(filterNames);

        if (settings.cluster.index < 0 || settings.cluster.index > filterNames.size())
        {
            cerr << "***Error: " << settings.cluster.index << " not a valid magnitude index.  Choose 0, 1,or 2.***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (1);
        }

        if (systems.size() > 1)
        {
            double logFieldStarLikelihood = 0.0;

            for (size_t filt = 0; filt < filterNames.size(); filt++)
            {
                logFieldStarLikelihood -= log (filterPriorMax.at(filt) - filterPriorMin.at(filt));
            }
            fsLike = exp (logFieldStarLikelihood);
        }
        else
        {
            fsLike = 0;
        }
    }

    // Begin initChain
    {
        for (auto system : systems)
        {
            system.clustStarProposalDens = system.clustStarPriorDens;   // Use prior prob of being clus star

            if (system.observedStatus == WD)
            {
                system.setMassRatio(0.0);
            }
        }
    }
    // end initChain

    cout << "Bayesian Analysis of Stellar Evolution" << endl;

    // Assuming fsLike doesn't change, this is the "global" logPost function
    auto logPostFunc = std::bind(&MpiMcmcApplication::logPostStep, this, _1, fsLike);

    // Run Burnin
    std::ofstream burninFile(ctrl.clusterFilename + ".burnin");
    if (!burninFile)
    {
        cerr << "***Error: File " << ctrl.clusterFilename + ".burnin" << " was not available for writing.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    printHeader (burninFile, ctrl.priorVar);

    Chain burninChain(static_cast<uint32_t>(std::uniform_int_distribution<>()(gen)), ctrl.priorVar, clust, burninFile);

    try
    {
        int adaptiveBurnIter = 0;
        bool acceptedOne = false;

        // Stage 1 burnin
        {
            cout << "\nRunning Stage 1 burnin..." << flush;

            auto proposalFunc = std::bind(&MpiMcmcApplication::propClustBigSteps, this, _1, std::cref(ctrl), std::cref(stepSize));
            burninChain.run(proposalFunc, logPostFunc, settings.mpiMcmc.adaptiveBigSteps);

            cout << " Complete (acceptanceRatio = " << burninChain.acceptanceRatio() << ")" << endl;
        }

        cout << "\nRunning Stage 2 (adaptive) burnin..." << endl;

        // Run adaptive burnin (stage 2)
        // -----------------------------
        // Exits after two consecutive trials with a 20% < x < 40% acceptanceRatio,
        // or after the numeber of iterations exceeds the maximum burnin iterations.
        do
        {
            // Convenience variable
            const int trialIter = settings.mpiMcmc.trialIter;

            // Increment the number of iterations we've gone through
            // If the step sizes aren't converging on an acceptable acceptance ratio, this can kick us out of the burnin
            adaptiveBurnIter += trialIter;

            // Reset the Chain for the next stage of the adaptive burnin
            burninChain.reset();

            std::function<Cluster(Cluster)> proposalFunc = std::bind(&MpiMcmcApplication::propClustIndep, this, _1, std::cref(ctrl), std::cref(stepSize), 5);

            if (settings.mpiMcmc.bigStepBurnin)
            {
                // Run big steps for the entire trial
                burninChain.run(proposalFunc, logPostFunc, trialIter);
            }
            else
            {
                // Run big steps for half the trial
                burninChain.run(proposalFunc, logPostFunc, trialIter / 2);

                // Then run smaller (but not the smallest) steps for the second half
                proposalFunc = std::bind(&MpiMcmcApplication::propClustIndep, this, _1, std::cref(ctrl), std::cref(stepSize), 1);
                burninChain.run(proposalFunc, logPostFunc, trialIter / 2);
            }

            double acceptanceRatio = burninChain.acceptanceRatio();

            if (acceptanceRatio <= 0.4 && acceptanceRatio >= 0.2)
            {
                if (acceptedOne)
                {
                    cout << "  Leaving adaptive burnin early with an acceptance ratio of " << acceptanceRatio << " (iteration " << adaptiveBurnIter << ")" << endl;
                    break;
                }
                else
                {
                    cout << "    Acceptance ratio: " << acceptanceRatio << ". Trying for trend." << endl;
                    acceptedOne = true;
                }
            }
            else
            {
                acceptedOne = false;
                cout << "    Acceptance ratio: " << acceptanceRatio << ". Retrying." << endl;

                scaleStepSizes(stepSize, burninChain.acceptanceRatio()); // Adjust step sizes
            }
        } while (adaptiveBurnIter < ctrl.burnIter);

        // Stage 3 burnin
        // Make sure and pull the covariance matrix before resetting the burninChain
        auto proposalFunc = std::bind(&MpiMcmcApplication::propClustCorrelated, this, _1, std::cref(ctrl), burninChain.makeCholeskyDecomp());

        burninChain.reset();

        cout << "\nRunning Stage 3 burnin... " << flush;

        burninChain.run(proposalFunc, logPostFunc, settings.mpiMcmc.stage3Iter);

        cout << " Complete (acceptanceRatio = " << burninChain.acceptanceRatio() << ")" << endl;

        burninFile.close();
    }
    catch (...)
    {
        burninFile.close();
        throw;
    }

    // Main run
    ofstream resultFile(ctrl.clusterFilename);
    if (!resultFile)
    {
        cerr << "***Error: File " << ctrl.clusterFilename << " was not available for writing.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    try
    {
        printHeader (resultFile, ctrl.priorVar);

        auto proposalFunc = std::bind(&MpiMcmcApplication::propClustCorrelated, this, _1, std::cref(ctrl), burninChain.makeCholeskyDecomp());

        Chain mainChain(static_cast<uint32_t>(std::uniform_int_distribution<>()(gen)), ctrl.priorVar, mainClust, resultFile);

        mainChain.run(proposalFunc, logPostFunc, ctrl.nIter, ctrl.thin);

        cout << "\nAcceptance ratio: " << mainChain.acceptanceRatio() << endl;

        resultFile.close();
    }
    catch (...)
    {
        resultFile.close();
        throw;
    }

    return 0;
}


void MpiMcmcApplication::scaleStepSizes (array<double, NPARAMS> &stepSize, double acceptanceRatio)
{
    function<double(double)> scaleFactor = [](double acceptanceRatio) {
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
    };

    for (int p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            stepSize.at(p) *= scaleFactor(acceptanceRatio);
        }
    }
}

Cluster MpiMcmcApplication::propClustBigSteps (Cluster clust, struct ifmrMcmcControl const &ctrl, array<double, NPARAMS> const &stepSize)
{
    return propClustIndep(clust, ctrl, stepSize, 25.0);
}

Cluster MpiMcmcApplication::propClustIndep (Cluster clust, struct ifmrMcmcControl const &ctrl, array<double, NPARAMS> const &stepSize, double scale)
{
    /* DOF defined in densities.h */
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            clust.setParam(p, clust.getParam(p) + sampleT (gen, scale * stepSize.at(p) * stepSize.at(p)));
        }
    }

    return clust;
}

Cluster MpiMcmcApplication::propClustCorrelated (Cluster clust, struct ifmrMcmcControl const &ctrl, Matrix<double, NPARAMS, NPARAMS> const &propMatrix)
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
                    corrProps.at(p) += propMatrix.at(p).at(k) * indepProps.at(k);
                }
            }
            clust.setParam(p, clust.getParam(p) + corrProps.at(p));
        }
    }

    return clust;
}

double MpiMcmcApplication::logPostStep(Cluster &propClust, double fsLike)
{
    mutex logPostMutex;
    double logPostProp;

    logPostProp = propClust.logPrior (evoModels);

    propClust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(propClust.feh, propClust.yyy, propClust.age);    // determine AGBt ZAMS mass, to find evol state

    /* loop over assigned stars */
    pool.parallelFor(systems.size(), [=,&logPostMutex,&logPostProp](int i)
    {
        double postClusterStar = 0.0;

        /* loop over all (mass1, mass ratio) pairs */
        if (systems.at(i).observedStatus == WD)
        {
            postClusterStar = 0.0;

            for (int j = 0; j < N_WD_MASS1; j++)
            {
                double tmpLogPost;
                StellarSystem wd(systems.at(i));
                wd.primary.mass = wdGridMass(j);
                wd.setMassRatio(0.0);

                try
                {
                    tmpLogPost = wd.logPost(propClust, evoModels);
                    tmpLogPost += log ((propClust.getM_wd_up() - MIN_MASS1) / (double) N_WD_MASS1);

                    postClusterStar +=  exp (tmpLogPost);
                }
                catch ( WDBoundsError &e )
                {
                    if (settings.verbose)
                        cerr << e.what() << endl;
                }
            }
        }
        else
        {
            try
            {
                /* marginalize over isochrone */
                postClusterStar = margEvolveWithBinary (propClust, systems.at(i), evoModels);
            }
            catch ( WDBoundsError &e )
            {
                if (settings.verbose)
                    cerr << e.what() << endl;
            }
        }

        postClusterStar *= systems.at(i).clustStarPriorDens;

        /* marginalize over field star status */
        std::lock_guard<mutex> lk(logPostMutex);
        logPostProp += log ((1.0 - systems.at(i).clustStarPriorDens) * fsLike + postClusterStar);
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
