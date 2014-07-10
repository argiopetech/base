#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/format.hpp>

#include "Chain.hpp"

#include "constants.hpp"
#include "densities.hpp"
#include "Matrix.hpp"
#include "Model.hpp"
#include "samplers.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "WhiteDwarf.hpp"
#include "Utility.hpp"

using std::array;
using std::string;
using std::istringstream;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::unique_ptr;

struct ifmrGridControl
{
    double initialAge;
    array<double, NPARAMS> priorMean, priorVar;
    double minMag;
    double maxMag;
    int iStart;
    int modelSet;
    vector<double> filterPriorMin, filterPriorMax;
    int nSamples;
    array<double, NPARAMS> start; /* starting points for grid evaluations */
    array<double, NPARAMS> end;   /* end points for grid evaluations */
};

/* For posterior evaluation on a grid */
struct clustPar
{
    clustPar(double age, double feh, double modulus, double absorption, double ifmrIntercept, double ifmrSlope, double ifmrQuadCoef)
        : age(age), FeH(feh), modulus(modulus), absorption(absorption), ifmrIntercept(ifmrIntercept), ifmrSlope(ifmrSlope), ifmrQuadCoef(ifmrQuadCoef)
    {}

    double age;
    double FeH;
    double modulus;
    double absorption;
    double ifmrIntercept;
    double ifmrSlope;
    double ifmrQuadCoef;
};

typedef struct
{
    vector<double> obsPhot, variance;
    double clustStarPriorDens;  /* cluster membership prior probability */
} obsStar;

/* Used in densities.c. */
array<double, NPARAMS> priorMean, priorVar;

/*
 * read control parameters from input stream
 */
static void initIfmrGridControl (Chain *mc, Model &evoModels, struct ifmrGridControl *ctrl, Settings &s)
{
    if (s.whiteDwarf.wdModel == WdModel::MONTGOMERY)
    {
        ctrl->priorMean.at(CARBONICITY) = mc->clust.carbonicity = mc->clust.priorMean.at(CARBONICITY) = s.cluster.carbonicity;
        ctrl->priorVar.at(CARBONICITY) = s.cluster.sigma.carbonicity;
    }
    else
    {
        ctrl->priorMean.at(CARBONICITY) = mc->clust.carbonicity = mc->clust.priorMean.at(CARBONICITY) = 0.0;
        ctrl->priorVar.at(CARBONICITY) = 0.0;
    }


    ctrl->priorMean.at(FEH) = s.cluster.Fe_H;
    ctrl->priorVar.at(FEH) = s.cluster.sigma.Fe_H;
    if (ctrl->priorVar.at(FEH) < 0.0)
    {
        ctrl->priorVar.at(FEH) = 0.0;
    }

    ctrl->priorMean.at(MOD) = s.cluster.distMod;
    ctrl->priorVar.at(MOD) = s.cluster.sigma.distMod;
    if (ctrl->priorVar.at(MOD) < 0.0)
    {
        ctrl->priorVar.at(MOD) = 0.0;
    }

    ctrl->priorMean.at(ABS) = s.cluster.Av;
    ctrl->priorVar.at(ABS) = s.cluster.sigma.Av;
    if (ctrl->priorVar.at(ABS) < 0.0)
    {
        ctrl->priorVar.at(ABS) = 0.0;
    }

    ctrl->initialAge = s.cluster.logClusAge;
    ctrl->priorVar.at(AGE) = 1.0;

    ctrl->priorVar.at(IFMR_INTERCEPT) = 1.0;
    ctrl->priorVar.at(IFMR_SLOPE) = 1.0;
    if (evoModels.IFMR >= 9)
        ctrl->priorVar.at(IFMR_QUADCOEF) = 1.0;
    else
        ctrl->priorVar.at(IFMR_QUADCOEF) = 0.0;

    // copy values to global variables
    priorVar.at(AGE) = ctrl->priorVar.at(AGE);
    priorVar.at(FEH) = ctrl->priorVar.at(FEH);
    priorVar.at(MOD) = ctrl->priorVar.at(MOD);
    priorVar.at(ABS) = ctrl->priorVar.at(ABS);
    priorVar.at(IFMR_INTERCEPT) = ctrl->priorVar.at(IFMR_INTERCEPT);
    priorVar.at(IFMR_SLOPE) = ctrl->priorVar.at(IFMR_SLOPE);
    priorVar.at(IFMR_QUADCOEF) = ctrl->priorVar.at(IFMR_QUADCOEF);

    priorMean.at(FEH) = ctrl->priorMean.at(FEH);
    priorMean.at(MOD) = ctrl->priorMean.at(MOD);
    priorMean.at(ABS) = ctrl->priorMean.at(ABS);

    /* prior values for linear IFMR */
    ctrl->priorMean.at(IFMR_SLOPE) = 0.08;
    ctrl->priorMean.at(IFMR_INTERCEPT) = 0.65;
    ctrl->priorMean.at(IFMR_QUADCOEF) = 0.0;
    priorMean.at(IFMR_SLOPE) = ctrl->priorMean.at(IFMR_SLOPE);
    priorMean.at(IFMR_INTERCEPT) = ctrl->priorMean.at(IFMR_INTERCEPT);
    priorMean.at(IFMR_QUADCOEF) = ctrl->priorMean.at(IFMR_QUADCOEF);

    /* open model file, choose model set, and load models */

    if (s.mainSequence.msRgbModel == MsModel::CHABHELIUM)
    {
        scanf ("%lf %lf", &ctrl->priorMean.at(YYY), &ctrl->priorVar.at(YYY));

        if (ctrl->priorVar.at(YYY) < 0.0)
        {
            ctrl->priorVar.at(YYY) = 0.0;
        }
    }
    else
    {
        ctrl->priorMean.at(YYY) = 0.0;
        ctrl->priorVar.at(YYY) = 0.0;
    }
    priorVar.at(YYY) = ctrl->priorVar.at(YYY);
    priorMean.at(YYY) = ctrl->priorMean.at(YYY);

    /* open files for reading (data) and writing */

    ctrl->minMag = s.cluster.minMag;
    ctrl->maxMag = s.cluster.maxMag;

    ctrl->iStart = 0;
} // initIfmrGridControl


/*
 * Read sampled params
 */
static void readSampledParams (struct ifmrGridControl *ctrl, vector<clustPar> &sampledPars, Model &evoModels, const Settings &s)
{
    string line;
    std::ifstream parsFile;
    parsFile.open(s.files.output + ".res");

    getline(parsFile, line); // Eat the header

    while (!parsFile.eof())
    {
        double newAge, newFeh, newMod, newAbs, newIInter, newISlope, newIQuad, ignore;
        newAge = newFeh = newMod = newAbs = newIInter = newISlope = newIQuad = 0.0;

        parsFile >> newAge
                 >> newFeh
                 >> newMod
                 >> newAbs;

        if (evoModels.IFMR >= 4)
        {
            parsFile >> newIInter
                     >> newISlope;
        }

        if (evoModels.IFMR >= 9)
        {
            parsFile >> newIQuad;
        }

        parsFile >> ignore; // logPost

        if (!parsFile.eof())
        {
            sampledPars.emplace_back(newAge, newFeh, newMod, newAbs, newIInter, newISlope, newIQuad);
        }
    }

    parsFile.close();

    ctrl->nSamples = sampledPars.size();
}


/*
 * Initialize chain
 */
static void initChain (Chain *mc, const struct ifmrGridControl *ctrl)
{
    int p;

    for (p = 0; p < NPARAMS; p++)
    {
        mc->acceptClust.at(p) = mc->rejectClust.at(p) = 0;
    }

    // If there is no beta in file, initialize everything to prior means
    mc->clust.feh = ctrl->priorMean.at(FEH);
    mc->clust.mod = ctrl->priorMean.at(MOD);
    mc->clust.abs = ctrl->priorMean.at(ABS);
    mc->clust.yyy = ctrl->priorMean.at(YYY);
    mc->clust.age = ctrl->initialAge;
    mc->clust.ifmrIntercept = ctrl->priorMean.at(IFMR_INTERCEPT);
    mc->clust.ifmrSlope = ctrl->priorMean.at(IFMR_SLOPE);
    mc->clust.ifmrQuadCoef = ctrl->priorMean.at(IFMR_QUADCOEF);
    mc->clust.mean.at(AGE) = ctrl->initialAge;
    mc->clust.mean.at(YYY) = ctrl->priorMean.at(YYY);
    mc->clust.mean.at(MOD) = ctrl->priorMean.at(MOD);
    mc->clust.mean.at(FEH) = ctrl->priorMean.at(FEH);
    mc->clust.mean.at(ABS) = ctrl->priorMean.at(ABS);
    mc->clust.mean.at(IFMR_INTERCEPT) = ctrl->priorMean.at(IFMR_INTERCEPT);
    mc->clust.mean.at(IFMR_SLOPE) = ctrl->priorMean.at(IFMR_SLOPE);
    mc->clust.mean.at(IFMR_QUADCOEF) = ctrl->priorMean.at(IFMR_QUADCOEF);

    for (auto star : mc->stars)
    {
        star.clustStarProposalDens = star.clustStarPriorDens;   // Use prior prob of being clus star

        star.primary.wdType = WdAtmosphere::DA;
        star.secondary.wdType = WdAtmosphere::DA;

        // find photometry for initial values of currentClust and mc->stars
        if (star.observedStatus == WD)
        {
            star.setMassRatio(0.0);
        }
    }
} // initChain


static bool acceptP (std::mt19937 &gen, const double current, const double proposed)
{
    if (std::isinf (proposed))
    {
//        cerr << "-Inf posterior proposed and rejected" << endl;
        return false;
    }

    double alpha = proposed - current;

    if (alpha >= 0) // Short circuit exit to the MH algorithm
    {
        return true;
    }

    double u = std::generate_canonical<double, 53>(gen);

    if (u < 1.e-15)
        u = 1.e-15;
    u = log (u);

    if (u < alpha)
    {
        return true;
    }
    else
    {
        return false;
    }
}


class Application
{
  private:
    Settings settings;
//    base::utility::ThreadPool pool;

  public:
    Application(Settings s)
        : settings(s)//, pool(s.threads)
    {}

    void run();
    std::tuple<double, double, double> sampleMass(const Cluster&, const Model&, const Isochrone&, std::mt19937&, StellarSystem, double, double);
};


std::tuple<double, double, double> Application::sampleMass(const Cluster &clust, const Model &evoModels, const Isochrone &isochrone, std::mt19937 &gen, StellarSystem star, const double massStepSize, const double massRatioStepSize)
{
    constexpr int iters = 12000;
    constexpr int burnIters = iters / 6;
    constexpr double scale = 10.0;
//     constexpr double massStepSize = 0.0025;
//     constexpr double massRatioStepSize = 0.053;
//     constexpr double massRatioStepSize = 0.006;

    double acceptedPosterior = -std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < iters; ++iter)
    {
        StellarSystem propStar = star;

        double tMass = propStar.primary.mass + sampleT (gen, (iter < burnIters ? scale : 1.0) * massStepSize * massStepSize);
        double tMassRatio = propStar.getMassRatio() + sampleT (gen, (iter < burnIters ? scale : 1.0) * massRatioStepSize * massRatioStepSize);

        if ((tMassRatio >= 0.0) && (tMassRatio <= 1.0) && (tMass > 0.1) && (tMass < clust.getM_wd_up()))
        {
            propStar.primary.mass = tMass;
            propStar.setMassRatio(tMassRatio);

            try
            {
                auto proposedPosterior = propStar.logPost(clust, evoModels, isochrone);

                if (acceptP(gen, acceptedPosterior, proposedPosterior))
                {
                    star = propStar;
                    acceptedPosterior = proposedPosterior;
                }
            }
            catch ( WDBoundsError &e )
            {
                // Go ahead and silence this error...
                cerr << e.what() << endl;
            }
        }
    }

    return std::make_tuple(star.primary.mass, star.getMassRatio(), acceptedPosterior);
}


void Application::run()
{
    Chain mc;
    struct ifmrGridControl ctrl;

    double fsLike;
    const double dMass1 = settings.sampleMass.deltaMass;
    const double dMassRatio = settings.sampleMass.deltaMassRatio;

    vector<clustPar> sampledPars;

    vector<std::pair<double, double>> masses;
    vector<double> memberships;

    std::mt19937 gen(settings.seed * uint32_t(2654435761)); // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)
    {
        const int warmupIter = 10000;

        cout << "Warming up generator..." << flush;

        for (auto i = 0; i < warmupIter; ++i)
        {
            std::generate_canonical<double, 10>(gen);
        }

        cout << " Done." << endl;
        cout << "Generated " << warmupIter << " values." << endl;
    }

    Model evoModels = makeModel(settings);

    initIfmrGridControl (&mc, evoModels, &ctrl, settings);

    {
        // Read photometry and calculate fsLike
        vector<double> filterPriorMin;
        vector<double> filterPriorMax;

        std::ifstream rData(settings.files.phot);

        if (!rData)
        {
            cerr << "***Error: Photometry file " << settings.files.phot << " was not found.***" << endl;
            cerr << ".at(Exiting...)" << endl;
            exit (1);
        }

        auto ret = base::utility::readPhotometry (rData, filterPriorMin, filterPriorMax, settings);
        auto filterNames = ret.first;
        mc.stars = ret.second;

        evoModels.restrictFilters(filterNames);

        if (settings.cluster.index < 0 || settings.cluster.index > filterNames.size())
        {
            cerr << "*** Error: " << settings.cluster.index << " not a valid magnitude index.  Choose 0 through " << filterNames.size() - 1 << " ***"<< endl;
            cerr << "(Exiting...)" << endl;
            exit (1);
        }

        double logFieldStarLikelihood = 0.0;

        for (size_t filt = 0; filt < filterNames.size(); filt++)
        {
            logFieldStarLikelihood -= log (filterPriorMax.at(filt) - filterPriorMin.at(filt));
        }

        fsLike = exp (logFieldStarLikelihood);
    }

    initChain (&mc, &ctrl);

    mc.clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

    readSampledParams (&ctrl, sampledPars, evoModels, settings);
    cout << "sampledPars.at(0).age = " << sampledPars.at(0).age << endl;

    /********** compile results *********/
    /*** now report sampled masses and parameters ***/

    masses.resize(mc.stars.size());
    memberships.resize(mc.stars.size());

    // Open the file
    string filename = settings.files.output + ".massSamples";

    std::ofstream massSampleFile;
    massSampleFile.open(filename);
    if (!massSampleFile)
    {
        cerr << "***Error: File " << filename << " was not available for writing.***" << endl;
        cerr << ".at(Exiting...)" << endl;
        exit (1);
    }

    filename += ".membership";

    std::ofstream membershipFile;
    membershipFile.open(filename);
    if (!membershipFile)
    {
        cerr << "***Error: File " << filename << " was not available for writing.***" << endl;
        cerr << ".at(Exiting...)" << endl;
        exit (1);
    }

    for (int m = 0; m < ctrl.nSamples; m++)
    {
        mc.clust.age = sampledPars.at(m).age;
        mc.clust.feh = sampledPars.at(m).FeH;
        mc.clust.mod = sampledPars.at(m).modulus;
        mc.clust.abs = sampledPars.at(m).absorption;

        if (evoModels.IFMR >= 4)
        {
            mc.clust.ifmrIntercept = sampledPars.at(m).ifmrIntercept;
            mc.clust.ifmrSlope = sampledPars.at(m).ifmrSlope;
        }

        if (evoModels.IFMR >= 9)
        {
            mc.clust.ifmrQuadCoef = sampledPars.at(m).ifmrQuadCoef;
        }

        /************ sample WD masses for different parameters ************/
        unique_ptr<Isochrone> isochrone(evoModels.mainSequenceEvol->deriveIsochrone(mc.clust.feh, mc.clust.yyy, mc.clust.age));

        auto starSize = mc.stars.size();

        for (decltype(starSize) i = 0; i < starSize; ++i)
        {
            auto sampleTuple = sampleMass(mc.clust, evoModels, *isochrone, gen, mc.stars.at(i), dMass1, dMassRatio);

            double postClusterStar = std::get<2>(sampleTuple);
            postClusterStar *= (mc.clust.getM_wd_up() - 0.15);

            masses.at(i) = std::pair<double, double>(std::get<0>(sampleTuple), std::get<1>(sampleTuple));

            memberships.at(i) = mc.stars.at(i).clustStarPriorDens * postClusterStar / (mc.stars.at(i).clustStarPriorDens * postClusterStar + (1.0 - mc.stars.at(i).clustStarPriorDens) * fsLike);
        }

        massSampleFile << boost::format("%10.6f") % sampledPars.at(m).age
                       << boost::format("%10.6f") % sampledPars.at(m).FeH
                       << boost::format("%10.6f") % sampledPars.at(m).modulus
                       << boost::format("%10.6f") % sampledPars.at(m).absorption;

        if (evoModels.IFMR >= 4)
        {
            massSampleFile << boost::format("%10.6f") % sampledPars.at(m).ifmrIntercept
                           << boost::format("%10.6f") % sampledPars.at(m).ifmrSlope;
        }

        if (evoModels.IFMR >= 9)
        {
            massSampleFile << boost::format("%10.6f") % sampledPars.at(m).ifmrQuadCoef;
        }

        for (auto mass : masses)
        {
            massSampleFile << boost::format("%10.6f") % mass.first
                           << boost::format("%10.6f") % mass.second;
        }

        for (auto membership : memberships)
        {
            membershipFile << boost::format("%10.6f") % membership;
        }

        massSampleFile << endl;
        membershipFile << endl;
    }

    massSampleFile.close();
    membershipFile.close();

    cout << "Completed successfully" << endl;
}

int main (int argc, char *argv[])
{
    Settings settings;

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

    Application(settings).run();

    return 0;
}
