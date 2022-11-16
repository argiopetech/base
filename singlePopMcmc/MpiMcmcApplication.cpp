#include <functional>
#include <iostream>
#include <memory>

#include "Chain.hpp"
#include "Cluster.hpp"
#include "IO/BackingStore.hpp"
#include "IO/Records.hpp"
#include "marg.hpp"
#include "mpiMcmc.hpp"
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
using std::unique_ptr;
using std::vector;

using namespace std::placeholders;

const double TWO_M_PI = 2 * M_PI;

void ensurePriors(const Settings &s, const Cluster &clust)
{
    // Carbonicity
    if (clust.carbonicity < 0.0)
        throw InvalidCluster("Low carbonicity");
    else if (clust.carbonicity > 1.0)
        throw InvalidCluster("High carbonicity");

    // Parallax
    // Bounded between 1 and 10^5 parsecs
    else if (s.modIsParallax && clust.mod < 0.00001)
        throw InvalidCluster("Low parallax");
    else if (s.modIsParallax && clust.mod > 1.0)
        throw InvalidCluster("High parallax");
}

MpiMcmcApplication::MpiMcmcApplication(Settings &s,
                                       SinglePopBackingStore *mcmcStore,
                                       FieldStarLikelihoodBackingStore *fieldStarLikelihood)
    : evoModels(makeModel(s)), settings(s)
    , gen(uint32_t(s.seed * uint32_t(2654435761))) // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)
    , mcmcStore(mcmcStore), fieldStarLikelihood(fieldStarLikelihood)
    , pool(s.threads)
{
    N_WD_MASS1 = s.wdInterpolationPower == 0 ? 8000 // Old default
                                             : s.wdInterpolationPower > 0 ? 64 << s.wdInterpolationPower // 2^(5 + wdInterpolationPower)
                                                                          : 64 >> -s.wdInterpolationPower;

    ctrl.priorVar.fill(0);

    clust.feh = settings.cluster.starting.Fe_H;
    clust.priorMean[FEH] = settings.cluster.priorMeans.Fe_H;
    ctrl.priorVar[FEH]   = settings.cluster.priorSigma.Fe_H;

    clust.mod = settings.cluster.starting.distMod;
    clust.priorMean[MOD] = settings.cluster.priorMeans.distMod;
    ctrl.priorVar[MOD]   = settings.cluster.priorSigma.distMod;

    clust.abs = settings.cluster.starting.Av;
    clust.priorMean[ABS] = fabs(settings.cluster.priorMeans.Av);
    ctrl.priorVar[ABS]   = settings.cluster.priorSigma.Av;

    clust.age = settings.cluster.starting.logAge;
    clust.priorMean[AGE] = settings.cluster.priorMeans.logAge;
    ctrl.priorVar[AGE]   = settings.cluster.priorSigma.logAge;

    clust.carbonicity = settings.cluster.starting.carbonicity;
    clust.priorMean[CARBONICITY] = settings.cluster.priorMeans.carbonicity;
    ctrl.priorVar[CARBONICITY]   = settings.cluster.priorSigma.carbonicity;

    clust.yyy = settings.cluster.starting.Y;
    clust.priorMean[YYY] = settings.cluster.priorMeans.Y;
    ctrl.priorVar[YYY]   = settings.cluster.priorSigma.Y;


    if (evoModels.IFMR <= 3 || evoModels.IFMR >= 12)
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

    clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

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
        ctrl.burnIter = settings.singlePopMcmc.burnIter;
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

    ctrl.nIter = settings.singlePopMcmc.maxIter;
    ctrl.thin = settings.singlePopMcmc.thin;

    std::copy(ctrl.priorVar.begin(), ctrl.priorVar.end(), clust.priorVar.begin());

    std::copy(s.singlePopMcmc.stepSize.begin(),
              s.singlePopMcmc.stepSize.end(),
              stepSize.begin());
}


// Initialize SSE memory for noBinaries
// THIS LEAKS MEMORY - FIXME
void MpiMcmcApplication::allocateSSEMem()
{
    using aligned_m128 = std::aligned_storage<16, 16>::type;

    if (! msSystems.empty())
    {
        // howManyFilts has to be a multiple of two to make the SSE code happy
        // Add 1 and round down.
        howManyFiltsAligned = ((msSystems.front().obsPhot.size() + 1) & ~0x1);
        howManyFilts        = msSystems.front().obsPhot.size();
        size_t howManyWeNeed       = msSystems.size() * howManyFiltsAligned;
        size_t howManyWeAlloc      = howManyWeNeed / 2;

        aligned_m128* talloc = new aligned_m128[howManyWeAlloc];
        sysVars = new(talloc) double[howManyWeNeed];

        talloc  = new aligned_m128[howManyWeAlloc];
        sysVar2 = new(talloc) double[howManyWeNeed];

        talloc = new aligned_m128[howManyWeAlloc];
        sysObs = new(talloc) double[howManyWeNeed];

    }

    int i = 0;

    for (auto s : msSystems)
    {
        for (size_t k = 0; k < howManyFiltsAligned; ++k, ++i)
        {
            if ((k < s.variance.size()) && (s.variance.at(k) > EPS))
            {
                sysVars[i] = s.variance.at(k) * clust.varScale;
                sysVar2[i] = __builtin_log (TWO_M_PI * sysVars[i]);

                // Removes a division from the noBinaries loop which turns a 14 cycle loop into a 3.8 cycle loop
                sysVars[i] = 1 / sysVars[i];

                sysObs [i] = s.obsPhot.at(k);
            }
            else
            {
                sysVars[i] = 0.0;
                sysVar2[i] = 0.0;

                sysObs [i] = 0.0;
            }
        }
    }
}

void MpiMcmcApplication::readPhotometry()
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

    for (auto r : ret.second)
    {
        if (r.observedStatus == StarStatus::MSRG)
        {
            // Everything goes into the main run
            msMainRun.push_back(r);

            if (r.useDuringBurnIn)
            {
                msSystems.push_back(r);
            }
        }
        else if (r.observedStatus == StarStatus::WD)
        {
            wdMainRun.push_back(r);

            if (r.useDuringBurnIn)
            {
                wdSystems.push_back(r);
            }
        }
        else
            cerr << "Found unsupported star in photometry, type '" << r.observedStatus << "'... Continuing anyway." << endl;
    }

    if ( msSystems.empty() && wdSystems.empty())
    {
        cerr << "No stars loaded... Exiting." << endl;
        exit(-1);
    }


    evoModels.restrictFilters(filterNames, settings.allowInvalidModels);

    if (settings.cluster.index < 0 || static_cast<size_t>(settings.cluster.index) > filterNames.size())
    {
        cerr << "***Error: " << settings.cluster.index << " not a valid magnitude index." << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    if ((msSystems.size() + wdSystems.size()) > 1)
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

void MpiMcmcApplication::verifyModelBounds()
{
    if (settings.verbose)
    {
        cout << std::setprecision(3) << std::fixed;
        cout << "Model boundaries are (" << evoModels.mainSequenceEvol->getMinAge()
             << ", " << evoModels.mainSequenceEvol->getMaxAge() << ") log years." << endl;

        cout << "Binaries are " << (settings.noBinaries ? "OFF" : "ON") << endl;
    }

    if (   settings.cluster.starting.logAge < evoModels.mainSequenceEvol->getMinAge()
        || settings.cluster.starting.logAge > evoModels.mainSequenceEvol->getMaxAge())
    {
        if (! settings.overrideBounds)
        {
            cerr << std::setprecision(3) << std::fixed
                 << "\n***Error: Starting value \"logAge\" (" << settings.cluster.starting.logAge
                 << ") is outside the model boundaries (" << evoModels.mainSequenceEvol->getMinAge()
                 << ", " << evoModels.mainSequenceEvol->getMaxAge()
                 << ").\n   Use the `--overrideBounds` flag if you really want this.\n[Exiting...]"
                 << endl;

            exit(-1);
        }
        else if (settings.verbose)
        {
            cout << std::setprecision(3) << std::fixed
                 << "\n***Warning: logAge (" << settings.cluster.starting.logAge
                 << ") is outside the model boundaries (" << evoModels.mainSequenceEvol->getMinAge()
                 << ", " << evoModels.mainSequenceEvol->getMaxAge()
                 << ").\n   Continuing due to `--overrideBounds` flag." << endl;
        }
    }
}

void MpiMcmcApplication::initChain()
{
    for (auto system : msSystems)
    {
        system.clustStarProposalDens = system.clustStarPriorDens;
    }

    for (auto system : wdSystems)
    {
        system.clustStarProposalDens = system.clustStarPriorDens;

        system.setMassRatio(0.0);
    }
}

void MpiMcmcApplication::stage1Burnin(Chain<Cluster>& chain, std::function<void(const Cluster&)>& checkPriors, std::function<double(Cluster&)>& logPostFunc)
{
    if ( settings.verbose )
        cout << "\nRunning Stage 1 burnin..." << flush;

    auto proposalFunc = std::bind(&MpiMcmcApplication::propClustBigSteps, this, _1, std::cref(ctrl), std::cref(stepSize));
    chain.run(AdaptiveMcmcStage::FixedBurnin, proposalFunc, logPostFunc, checkPriors, settings.singlePopMcmc.adaptiveBigSteps);

    if ( settings.verbose )
        cout << " Complete (acceptanceRatio = " << chain.acceptanceRatio() << ")" << endl;

    chain.reset(); // Reset the chain to forget this part of the burnin.
}

void MpiMcmcApplication::stage2Burnin(Chain<Cluster>& chain, std::function<void(const Cluster&)>& checkPriors, std::function<double(Cluster&)>& logPostFunc)
{
    if ( settings.verbose )
        cout << "\nRunning Stage 2 (adaptive) burnin..." << endl;

    int  adaptiveBurnIter = 0;
    bool acceptedOne = false;  // Keeps track of whether we have two accepted trials in a row
    bool acceptedOnce = false; // Keeps track of whether we have ever accepted a trial

    // Run adaptive burnin (stage 2)
    // -----------------------------
    // Exits after two consecutive trials with a 20% < x < 40% acceptanceRatio,
    // or after the numeber of iterations exceeds the maximum burnin iterations.
    do
    {
        // Convenience variable
        const int trialIter = settings.singlePopMcmc.trialIter;

        // Increment the number of iterations we've gone through
        // If the step sizes aren't converging on an acceptable acceptance ratio, this can kick us out of the burnin
        adaptiveBurnIter += trialIter;

        // Reset the ratio to determine step scaling for this trial
        chain.resetRatio();

        std::function<Cluster(Cluster)> proposalFunc = std::bind(&MpiMcmcApplication::propClustIndep, this, _1, std::cref(ctrl), std::cref(stepSize), 5);

        if (settings.singlePopMcmc.bigStepBurnin)
        {
            // Run big steps for the entire trial
            chain.run(AdaptiveMcmcStage::AdaptiveBurnin, proposalFunc, logPostFunc, checkPriors, trialIter);
        }
        else if (acceptedOnce)
        {
            proposalFunc = std::bind(&MpiMcmcApplication::propClustIndep, this, _1, std::cref(ctrl), std::cref(stepSize), 1);
            chain.run(AdaptiveMcmcStage::AdaptiveBurnin, proposalFunc, logPostFunc, checkPriors, trialIter / 2);
        }
        else
        {
            // Run big steps for half the trial
            chain.run(AdaptiveMcmcStage::AdaptiveBurnin, proposalFunc, logPostFunc, checkPriors, trialIter / 2);

            // Then run smaller steps for the second half
            proposalFunc = std::bind(&MpiMcmcApplication::propClustIndep, this, _1, std::cref(ctrl), std::cref(stepSize), 1);
            chain.run(AdaptiveMcmcStage::AdaptiveBurnin, proposalFunc, logPostFunc, checkPriors, trialIter / 2);
        }

        double acceptanceRatio = chain.acceptanceRatio();

        if (acceptanceRatio <= 0.4 && acceptanceRatio >= 0.2)
        {
            acceptedOnce = true;

            if (acceptedOne)
            {
                if ( settings.verbose )
                    cout << "  Leaving adaptive burnin early with an acceptance ratio of " << acceptanceRatio << " (iteration " << adaptiveBurnIter + settings.singlePopMcmc.adaptiveBigSteps << ")" << endl;

                break;
            }
            else
            {
                if ( settings.verbose )
                    cout << "    Acceptance ratio: " << acceptanceRatio << ". Trying for trend." << endl;
                acceptedOne = true;
            }
        }
        else
        {
            acceptedOne = false;

            if ( settings.verbose )
                cout << "    Acceptance ratio: " << acceptanceRatio << ". Retrying." << endl;

            scaleStepSizes(stepSize, chain.acceptanceRatio()); // Adjust step sizes
        }
    } while (adaptiveBurnIter < ctrl.burnIter);
}


void MpiMcmcApplication::stage3Burnin(Chain<Cluster>& chain, std::function<void(const Cluster&)>& checkPriors, std::function<double(Cluster&)>& logPostFunc)
{
    if ( settings.verbose )
        cout << "\nStarting adaptive run... " << flush;

    // Make sure and pull the covariance matrix before resetting the chain
    auto proposalFunc = std::bind(&MpiMcmcApplication::propClustCorrelated, this, _1, std::cref(ctrl), chain.makeCholeskyDecomp());

    chain.reset();

    chain.run(AdaptiveMcmcStage::AdaptiveMainRun, proposalFunc, logPostFunc, checkPriors, settings.singlePopMcmc.stage3Iter);

    if ( settings.verbose )
        cout << " Preliminary acceptanceRatio = " << chain.acceptanceRatio() << endl;
}


void MpiMcmcApplication::mainRun(Chain<Cluster>& chain, std::function<void(const Cluster&)>& checkPriors, std::function<double(Cluster&)>& logPostFunc)
{
    // Begin main run
    // Main run proceeds in increments of 1, adapting the covariance matrix after every increment
    for (auto iters = 0; iters < ctrl.nIter; ++iters)
    {
        auto proposalFunc = std::bind(&MpiMcmcApplication::propClustCorrelated, this, _1, std::cref(ctrl), chain.makeCholeskyDecomp());

        chain.run(AdaptiveMcmcStage::MainRun, proposalFunc, logPostFunc, checkPriors, 1, ctrl.thin);
    }
}


void MpiMcmcApplication::loadPreviousBurnin(string filename, Chain<Cluster>& chain)
{
    auto acceptStage = [](int stage) { return stage == (int)AdaptiveMcmcStage::AdaptiveBurnin; };
    auto samples = base::utility::readSampledParams(filename, evoModels, settings, acceptStage);

    for (auto sample : samples)
    {
        clust.age = sample.age;
        clust.yyy = sample.y;
        clust.feh = sample.feh;
        clust.mod = sample.distMod;
        clust.abs = sample.abs;
        clust.carbonicity = sample.carbonicity;
        clust.ifmrIntercept = sample.ifmrIntercept;
        clust.ifmrSlope = sample.ifmrSlope;
        clust.ifmrQuadCoef = sample.ifmrQuadCoef;

        auto logPost = sample.logPost;

        chain.storePars(sample.stage, clust, logPost);
    }
}


int MpiMcmcApplication::run()
{
    Chain<Cluster> chain(static_cast<uint32_t>(std::uniform_int_distribution<>()(gen)),
                         ctrl.priorVar, clust, *mcmcStore, settings.modIsParallax);


    verifyModelBounds();

    readPhotometry();

    fieldStarLikelihood->save({fsLike});

    // Burnin
    allocateSSEMem();
    initChain();

    std::function<double(Cluster&)> logPostFunc = std::bind(&MpiMcmcApplication::logPostStep, this, _1);
    std::function<void(const Cluster&)> checkPriors = std::bind(&ensurePriors, std::cref(settings), _1);

    if (settings.startWithBurnin.empty())
    {
        stage1Burnin(chain, checkPriors, logPostFunc);
        stage2Burnin(chain, checkPriors, logPostFunc);
    }
    else
    {
        loadPreviousBurnin(settings.startWithBurnin, chain);
    }

    if (!settings.stopAfterBurnin)
    {
        stage3Burnin(chain, checkPriors, logPostFunc);

        // Main run
        // Overwrite the burnin photometry set with the full photometry set
        // SSE memory must be reallocated to fit the potentially larger number of stars
        msSystems = msMainRun;
        wdSystems = wdMainRun;

        allocateSSEMem();

        mainRun(chain, checkPriors, logPostFunc);


        // Post-run
        if ( settings.verbose )
        {
            cout << "\nFinal acceptance ratio: " << chain.acceptanceRatio() << "\n" << endl;

            MsBoundsError::countSummary();
        }
    }
    else
        cout << "Ended after burnin due to `--stopAfterBurnin`" << endl;

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
    array<double, NPARAMS> tDraws;

    for (auto &d : tDraws)
    {
        d = sampleT (gen, 1.0);
    }

    for (int p = 0; p < NPARAMS; ++p)
    {
        if (ctrl.priorVar.at(p) > EPSILON)
        {
            double corrProps = 0;

            for (int k = 0; k <= p; ++k)
            {                           /* propMatrix is lower diagonal */
                if (ctrl.priorVar.at(k) > EPSILON)
                {
                    corrProps += propMatrix.at(p).at(k) * tDraws[k];
                }
            }

            clust.setParam(p, clust.getParam(p) + corrProps);
        }
    }

    return clust;
}

double MpiMcmcApplication::logPostStep(Cluster &propClust)
{
    double logPostProp = propClust.logPrior (evoModels);

    unique_ptr<Isochrone> isochrone(evoModels.mainSequenceEvol->deriveIsochrone(propClust.feh, propClust.yyy, propClust.age));

    // Handle WDs specially (this needs to be moved to a new function)
    if ( ! wdSystems.empty() )
    {
        Star s;
        vector<double> primaryMags, secondaryMags, daCombinedMags, dbCombinedMags;

        s.mass = 0.0;
        secondaryMags = s.getMags(propClust, evoModels, *isochrone);

        vector<double> post(wdSystems.size(), 0.0);

        const auto m_wd_up = propClust.getM_wd_up();
        const auto wdSize  = wdSystems.size();
        const auto wdMassPrior = __builtin_log ((m_wd_up - MIN_MASS1) / (double) N_WD_MASS1);

        for (size_t j = 0; j < N_WD_MASS1; ++j)
        {
            const double primaryMass = wdGridMass(j);
            const double logPrior    = propClust.logPriorMass (primaryMass);

            if (!settings.onlyWDs || primaryMass > isochrone->agbTipMass())
            {
                s.mass       = primaryMass;

                s.wdType     = WdAtmosphere::DA;
                primaryMags  = s.getMags(propClust, evoModels, *isochrone);
                daCombinedMags = StellarSystem::deriveCombinedMags(propClust, evoModels, *isochrone, primaryMags, secondaryMags, settings.modIsParallax);

                s.wdType     = WdAtmosphere::DB;
                primaryMags  = s.getMags(propClust, evoModels, *isochrone);
                dbCombinedMags = StellarSystem::deriveCombinedMags(propClust, evoModels, *isochrone, primaryMags, secondaryMags, settings.modIsParallax);

                for (size_t i = 0; i < wdSize; ++i)
                {
                    double tmpLogPost;

                    if (wdSystems[i].primary.wdType == WdAtmosphere::DA)
                    {
                        tmpLogPost  = wdSystems[i].logPost(propClust, evoModels, *isochrone, logPrior, daCombinedMags);
                    }
                    else
                    {
                        tmpLogPost  = wdSystems[i].logPost(propClust, evoModels, *isochrone, logPrior, dbCombinedMags);
                    }
                    tmpLogPost += wdMassPrior;

                    post[i] += __builtin_exp(tmpLogPost);
                }
            }
        }

        for (size_t i = 0; i < wdSize; ++i)
        {
            if (post[i] > 0.0)
            {
                post[i] *= wdSystems[i].clustStarPriorDens;
            }
            else
            {
                post[i] = 0.0;
            }

            /* marginalize over field star status */
            logPostProp += __builtin_log ((1.0 - wdSystems[i].clustStarPriorDens) * fsLike + post[i]);
        }
    }

    if ( ! msSystems.empty() )
    {
        auto msSize = msSystems.size();

        vector<double> post;

        if (settings.noBinaries)
            post = margEvolveNoBinaries (propClust, evoModels, *isochrone, pool, sysVars, sysVar2, sysObs, msSize, howManyFiltsAligned, howManyFilts, settings.modIsParallax, settings.eepInterpolationPower);
        else
            post = margEvolveWithBinary (propClust, msSystems, evoModels, *isochrone, pool, settings.modIsParallax, settings.eepInterpolationPower);

        for (size_t i = 0; i < msSize; ++i)
        {
            if (post[i] > 0.0)
            {
                post[i] *= msSystems[i].clustStarPriorDens;
            }
            else
            {
                post[i] = 0.0;
            }

            // marginalize over field star status
            logPostProp += __builtin_log ((1.0 - msSystems[i].clustStarPriorDens) * fsLike + post[i]);
        }
    }

    return logPostProp;
}

double MpiMcmcApplication::wdGridMass (int point) const
{
    // I think this only gets calculated once, but I can't quite be sure...
    // It may not even matter, but every little bit helps, eh? And hey, preoptimization...
    static const double dWdMass1 = (settings.whiteDwarf.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1;

    return MIN_MASS1 + (dWdMass1 * point);
}
