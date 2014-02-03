#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/format.hpp>

#include "Chain.hpp"

#include "constants.hpp"
#include "densities.hpp"
#include "evolve.hpp"
#include "Matrix.hpp"
#include "Model.hpp"
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

struct ifmrGridControl
{
    double initialAge;
    double priorMean[NPARAMS];
    double priorVar[NPARAMS];
    double minMag;
    double maxMag;
    int iMag;
    int iStart;
    int modelSet;
    double filterPriorMin[FILTS];
    double filterPriorMax[FILTS];
    int numFilts;
    int nSamples;
    double start[NPARAMS];      /* starting points for grid evaluations */
    double end[NPARAMS];        /* end points for grid evaluations */
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
    double obsPhot[FILTS];
    double variance[FILTS];
    double clustStarPriorDens;  /* cluster membership prior probability */
} obsStar;

/* declare global variables */
array<double, FILTS> filterPriorMin;
array<double, FILTS> filterPriorMax;

/* Used in densities.c. */
double priorMean[NPARAMS], priorVar[NPARAMS];

/* Used by a bunch of different functions. */
vector<int> filters;

/*
 * read control parameters from input stream
 */
static void initIfmrGridControl (Chain *mc, Model &evoModels, struct ifmrGridControl *ctrl, Settings &s)
{
    ctrl->numFilts = 0;

    if (s.whiteDwarf.wdModel == WdModel::MONTGOMERY)
    {
        ctrl->priorMean[CARBONICITY] = mc->clust.carbonicity = mc->clust.priorMean[CARBONICITY] = s.cluster.carbonicity;
        ctrl->priorVar[CARBONICITY] = s.cluster.sigma.carbonicity;
    }
    else
    {
        ctrl->priorMean[CARBONICITY] = mc->clust.carbonicity = mc->clust.priorMean[CARBONICITY] = 0.0;
        ctrl->priorVar[CARBONICITY] = 0.0;
    }


    ctrl->priorMean[FEH] = s.cluster.Fe_H;
    ctrl->priorVar[FEH] = s.cluster.sigma.Fe_H;
    if (ctrl->priorVar[FEH] < 0.0)
    {
        ctrl->priorVar[FEH] = 0.0;
    }

    ctrl->priorMean[MOD] = s.cluster.distMod;
    ctrl->priorVar[MOD] = s.cluster.sigma.distMod;
    if (ctrl->priorVar[MOD] < 0.0)
    {
        ctrl->priorVar[MOD] = 0.0;
    }

    ctrl->priorMean[ABS] = s.cluster.Av;
    ctrl->priorVar[ABS] = s.cluster.sigma.Av;
    if (ctrl->priorVar[ABS] < 0.0)
    {
        ctrl->priorVar[ABS] = 0.0;
    }

    ctrl->initialAge = s.cluster.logClusAge;
    ctrl->priorVar[AGE] = 1.0;

    ctrl->priorVar[IFMR_INTERCEPT] = 1.0;
    ctrl->priorVar[IFMR_SLOPE] = 1.0;
    if (evoModels.IFMR >= 9)
        ctrl->priorVar[IFMR_QUADCOEF] = 1.0;
    else
        ctrl->priorVar[IFMR_QUADCOEF] = 0.0;

    // copy values to global variables
    priorVar[AGE] = ctrl->priorVar[AGE];
    priorVar[FEH] = ctrl->priorVar[FEH];
    priorVar[MOD] = ctrl->priorVar[MOD];
    priorVar[ABS] = ctrl->priorVar[ABS];
    priorVar[IFMR_INTERCEPT] = ctrl->priorVar[IFMR_INTERCEPT];
    priorVar[IFMR_SLOPE] = ctrl->priorVar[IFMR_SLOPE];
    priorVar[IFMR_QUADCOEF] = ctrl->priorVar[IFMR_QUADCOEF];

    priorMean[FEH] = ctrl->priorMean[FEH];
    priorMean[MOD] = ctrl->priorMean[MOD];
    priorMean[ABS] = ctrl->priorMean[ABS];

    /* prior values for linear IFMR */
    ctrl->priorMean[IFMR_SLOPE] = 0.08;
    ctrl->priorMean[IFMR_INTERCEPT] = 0.65;
    ctrl->priorMean[IFMR_QUADCOEF] = 0.0;
    priorMean[IFMR_SLOPE] = ctrl->priorMean[IFMR_SLOPE];
    priorMean[IFMR_INTERCEPT] = ctrl->priorMean[IFMR_INTERCEPT];
    priorMean[IFMR_QUADCOEF] = ctrl->priorMean[IFMR_QUADCOEF];

    /* open model file, choose model set, and load models */

    if (s.mainSequence.msRgbModel == MsModel::CHABHELIUM)
    {
        scanf ("%lf %lf", &ctrl->priorMean[YYY], &ctrl->priorVar[YYY]);

        if (ctrl->priorVar[YYY] < 0.0)
        {
            ctrl->priorVar[YYY] = 0.0;
        }
    }
    else
    {
        ctrl->priorMean[YYY] = 0.0;
        ctrl->priorVar[YYY] = 0.0;
    }
    priorVar[YYY] = ctrl->priorVar[YYY];
    priorMean[YYY] = ctrl->priorMean[YYY];

    /* open files for reading (data) and writing */

    ctrl->minMag = s.cluster.minMag;
    ctrl->maxMag = s.cluster.maxMag;
    ctrl->iMag = s.cluster.index;
    if (ctrl->iMag < 0 || ctrl->iMag > FILTS)
    {
        cerr << "***Error: " << ctrl->iMag << " not a valid magnitude index.  Choose 0, 1,or 2.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    ctrl->iStart = 0;

    /* Initialize filter prior mins and maxes */
    int j;

    for (j = 0; j < FILTS; j++)
    {
        ctrl->filterPriorMin[j] = 1000;
        ctrl->filterPriorMax[j] = -1000;
    }
} // initIfmrGridControl

void readCmdData (Chain &mc, struct ifmrGridControl &ctrl, const Model &evoModels, const Settings &s)
{
    string line, pch;

    std::ifstream rData;
    rData.open(s.files.phot);
    if (!rData)
    {
        cerr << "***Error: Photometry file " << s.files.phot << " was not found.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    //Parse the header of the file to determine which filters are being used
    getline(rData, line);  // Read in the header line

    istringstream header(line); // Ignore the first token (which is "id")

    header >> pch;

    while (!header.eof())
    {
        header >> pch;

        if (pch == "sig")
            break;                      // and check to see if they are 'sig'.  If they are, there are no more filters

        for (int filt = 0; filt < FILTS; filt++)
        {                               // Otherwise check to see what this filter's name is
            if (pch == evoModels.filterSet->getFilterName(filt))
            {
                ctrl.numFilts++;
                filters.push_back(filt);
                const_cast<Model&>(evoModels).numFilts++;
                break;
            }
        }
    }

    // This loop reads in photometry data
    // It also reads a best guess for the mass
    mc.stars.clear();

    while (!rData.eof())
    {
        getline(rData, line);

        if (rData.eof())
            break;

        mc.stars.push_back(Star(line, ctrl.numFilts));

        for (int i = 0; i < ctrl.numFilts; i++)
        {
            if (mc.stars.back().obsPhot[i] < ctrl.filterPriorMin[i])
            {
                ctrl.filterPriorMin[i] = mc.stars.back().obsPhot[i];
            }

            if (mc.stars.back().obsPhot[i] > ctrl.filterPriorMax[i])
            {
                ctrl.filterPriorMax[i] = mc.stars.back().obsPhot[i];
            }

            filterPriorMin[i] = ctrl.filterPriorMin[i];
            filterPriorMax[i] = ctrl.filterPriorMax[i];
        }

        if (!(mc.stars.back().status[0] == 3 || (mc.stars.back().obsPhot[ctrl.iMag] >= ctrl.minMag && mc.stars.back().obsPhot[ctrl.iMag] <= ctrl.maxMag)))
        {
            mc.stars.pop_back();
        }
    }

    rData.close();
} /* readCmdData */

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
        mc->acceptClust[p] = mc->rejectClust[p] = 0;
    }

    // If there is no beta in file, initialize everything to prior means
    mc->clust.feh = ctrl->priorMean[FEH];
    mc->clust.mod = ctrl->priorMean[MOD];
    mc->clust.abs = ctrl->priorMean[ABS];
    mc->clust.yyy = ctrl->priorMean[YYY];
    mc->clust.age = ctrl->initialAge;
    mc->clust.ifmrIntercept = ctrl->priorMean[IFMR_INTERCEPT];
    mc->clust.ifmrSlope = ctrl->priorMean[IFMR_SLOPE];
    mc->clust.ifmrQuadCoef = ctrl->priorMean[IFMR_QUADCOEF];
    mc->clust.mean[AGE] = ctrl->initialAge;
    mc->clust.mean[YYY] = ctrl->priorMean[YYY];
    mc->clust.mean[MOD] = ctrl->priorMean[MOD];
    mc->clust.mean[FEH] = ctrl->priorMean[FEH];
    mc->clust.mean[ABS] = ctrl->priorMean[ABS];
    mc->clust.mean[IFMR_INTERCEPT] = ctrl->priorMean[IFMR_INTERCEPT];
    mc->clust.mean[IFMR_SLOPE] = ctrl->priorMean[IFMR_SLOPE];
    mc->clust.mean[IFMR_QUADCOEF] = ctrl->priorMean[IFMR_QUADCOEF];

    int i;

    for (auto star : mc->stars)
    {
        star.meanMassRatio = 0.0;
        star.isFieldStar = 0;
        star.clustStarProposalDens = star.clustStarPriorDens;   // Use prior prob of being clus star
        star.UStepSize = 0.001; // within factor of ~2 for most main sequence stars
        star.massRatioStepSize = 0.001;

        for (i = 0; i < NPARAMS; i++)
        {
            star.beta[i][0] = 0.0;
            star.beta[i][1] = 0.0;
        }
        star.betaMassRatio[0] = 0.0;
        star.betaMassRatio[1] = 0.0;
        star.meanU = 0.0;
        star.varU = 0.0;

        for (i = 0; i < 2; i++)
            star.wdType[i] = WdAtmosphere::DA;

        for (i = 0; i < ctrl->numFilts; i++)
        {
            star.photometry[i] = 0.0;
        }

        // find photometry for initial values of currentClust and mc->stars
        if (star.status[0] == WD)
        {
            star.UStepSize = 0.05;      // use larger initial step size for white dwarfs
            star.massRatio = 0.0;
        }
    }
} // initChain


class Application
{
  private:
    Settings settings;
    base::utility::ThreadPool pool;

  public:
    Application(Settings s)
        : settings(s), pool(s.threads)
    {}

    void run();
    std::tuple<double, double, double> sampleMass(const double, const Cluster&, const Model&, const double, const double, const double, const double, Star);
};

// O(3n²)
std::tuple<double, double, double> Application::sampleMass(const double U, const Cluster &clust, const Model &evoModels, const double baseMass, const double maxMass, const double deltaPrimaryMass, const double deltaMassRatio, Star star)
{
    double maxLogPost = std::numeric_limits<double>::min();
    double cumulative = 0.0, denom = 0.0, postClusterStar = 0.0;
    int primaryMassIndex = 0, massRatioIndex = 0;

    const int maxPrimaryIndex = std::ceil((maxMass - baseMass) / deltaPrimaryMass);
    const int maxRatioIndex   = std::ceil(1.0 / deltaMassRatio);

    Vatrix<double> logPosts;
    logPosts.reserve( maxPrimaryIndex + 1 );

    // O(n²)
    for (int i = 0; i <= maxPrimaryIndex; ++i)
    {
        const double primaryMass = baseMass + (deltaPrimaryMass * i);

        logPosts.emplace_back();
        logPosts.back().reserve( maxRatioIndex + 1 );
        
        for (int j = 0; j <= maxRatioIndex; ++j)
        {
            const double massRatio = (deltaMassRatio * j);

            star.U = primaryMass;
            star.massRatio = massRatio;

            try
            {
                array<double, 2> ltau;
                array<double, FILTS> globalMags;
                evolve (clust, evoModels, globalMags, filters, star, ltau);

                auto thisLogPost = logPost1Star(star, clust, evoModels, filterPriorMin, filterPriorMax);
                logPosts.back().push_back( thisLogPost );

                postClusterStar += exp(thisLogPost);

                maxLogPost = std::max(maxLogPost, thisLogPost);
            }
            catch ( WDBoundsError &e )
            {
                // Go ahead and silence this error...
                cerr << e.what() << endl;
    
                logPosts.back().push_back(std::numeric_limits<double>::min());
            }
        }
    }

    // O(n²)
    for (auto v : logPosts)
    {
        for (auto logPost : v)
        {
            denom += exp(logPost - maxLogPost);
        }
    }

    // O(n²)
    while (cumulative <= U)
    {
        cumulative += exp(logPosts.at(primaryMassIndex).at(massRatioIndex) - maxLogPost) / denom;
        ++massRatioIndex;
        
        if (massRatioIndex > maxRatioIndex)
        {
            massRatioIndex = 0;
            ++primaryMassIndex; // This should never exceed primaryMasses.size() before cumulative > U.
        }
    }

    return std::make_tuple(baseMass + (primaryMassIndex * deltaPrimaryMass), massRatioIndex * deltaMassRatio, postClusterStar);
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

    readCmdData (mc, ctrl, evoModels, settings);

    evoModels.numFilts = ctrl.numFilts;

    initChain (&mc, &ctrl);

    mc.clust.M_wd_up = settings.whiteDwarf.M_wd_up;

    for (decltype(mc.stars.size()) i = 0; i < mc.stars.size(); i++)
    {
        mc.stars.at(i).isFieldStar = 0;
    }

    double logFieldStarLikelihood = 0.0;

    for (int filt = 0; filt < ctrl.numFilts; filt++)
    {
        logFieldStarLikelihood -= log (ctrl.filterPriorMax[filt] - ctrl.filterPriorMin[filt]);
    }
    fsLike = exp (logFieldStarLikelihood);

    readSampledParams (&ctrl, sampledPars, evoModels, settings);
    cout << "sampledPars[0].age = " << sampledPars.at(0).age << endl;

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
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    filename += ".membership";

    std::ofstream membershipFile;
    membershipFile.open(filename);
    if (!membershipFile)
    {
        cerr << "***Error: File " << filename << " was not available for writing.***" << endl;
        cerr << "[Exiting...]" << endl;
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
       mc.clust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, mc.clust.feh, mc.clust.yyy, mc.clust.age);

        for (int i = 0; i < mc.stars.size(); ++i)
        {
            auto sampleTuple = sampleMass(std::generate_canonical<double, 53>(gen), mc.clust, evoModels, 0.15, mc.clust.M_wd_up, dMass1, dMassRatio, mc.stars.at(i));

            double postClusterStar = std::get<2>(sampleTuple);
            postClusterStar *= (mc.clust.M_wd_up - 0.15);

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
