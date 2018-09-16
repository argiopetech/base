#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

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
using std::stringstream;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::unique_ptr;

using namespace std::placeholders;

/* For posterior evaluation on a grid */
struct clustPar
{
    clustPar(double age, double y, double feh, double modulus, double absorption, double carbonicity, double ifmrIntercept, double ifmrSlope, double ifmrQuadCoef, double logPost)
        : age(age), y(y), feh(feh), distMod(modulus), abs(absorption), carbonicity(carbonicity), ifmrIntercept(ifmrIntercept), ifmrSlope(ifmrSlope), ifmrQuadCoef(ifmrQuadCoef), logPost(logPost)
    {}

    double age;
    double y;
    double feh;
    double distMod;
    double abs;
    double carbonicity;
    double ifmrIntercept;
    double ifmrSlope;
    double ifmrQuadCoef;
    double logPost;
};


/*
 * Read sampled params
 */
static vector<clustPar> readSampledParams (Model &evoModels, const Settings &s)
{
    vector<clustPar> sampledPars;

    if (s.run <= 0)
    {
        string line;

        std::ifstream parsFile;
        parsFile.open(s.files.output + ".res");

        bool hasY, hasCarbonicity;

        getline(parsFile, line); // Parse header

        {
            string sin;
            stringstream in(line);

            in >> sin  // logAge
               >> sin; // Y?

            if (sin == "Y")
            {
                hasY = true;

                in >> sin  // FeH
                   >> sin  // Abs
                   >> sin  // DistMod
                   >> sin; // Carbonicity?

                if (sin == "carbonicity")
                    hasCarbonicity = true;
                else
                    hasCarbonicity = false;
            }
            else
            {
                hasY = false;

                // This one skips FeH (because it was already read instead of Y)
                in >> sin  // Abs
                   >> sin  // DistMod
                   >> sin; // Carbonicity?

                if (sin == "carbonicity")
                    hasCarbonicity = true;
                else
                    hasCarbonicity = false;
            }
        }

        while (getline(parsFile, line))
        {
            stringstream in(line);

            double newAge, newY, newFeh, newMod, newAbs, newCarbonicity, newIInter, newISlope, newIQuad, newLogPost;

            in >> newAge;

            if (hasY)
                in >> newY;
            else
                newY = s.cluster.starting.Y;

            in >> newFeh
               >> newMod
               >> newAbs;

            if (hasCarbonicity)
                in >> newCarbonicity;
            else
                newCarbonicity = s.cluster.starting.carbonicity;

            if (evoModels.IFMR >= 4)
            {
                in >> newIInter
                   >> newISlope;
            }

            if (evoModels.IFMR >= 9)
            {
                in >> newIQuad;
            }

            in >> newLogPost;

            sampledPars.emplace_back(newAge, newY, newFeh, newMod, newAbs, newCarbonicity, newIInter, newISlope, newIQuad, newLogPost);
        }

        parsFile.close();
    }
    else
    {
        auto records = base::utility::readSinglePopMainRunResFromDB(s.run, s.files.output);

        for (auto r : records)
        {
            sampledPars.emplace_back(r.clust.age, r.clust.yyy, r.clust.feh, r.clust.mod,
                                     r.clust.abs, r.clust.carbonicity, r.clust.ifmrIntercept,
                                     r.clust.ifmrSlope, r.clust.ifmrQuadCoef, r.logPost);
        }
    }

    return sampledPars;
}


class Application
{
  private:
    Settings settings;

    std::mt19937 gen;

    Cluster clust;
    Model evoModels;

    const double massStepSize;
    const double massRatioStepSize;

  public:
    Application(Settings s)
        // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)
        : settings(s), gen(s.seed * uint32_t(2654435761)), evoModels(makeModel(s))
        , massStepSize(settings.sampleMass.deltaMass)
        , massRatioStepSize( s.noBinaries ? 0.0 : settings.sampleMass.deltaMassRatio )
    {
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

        clust.feh = settings.cluster.starting.Fe_H;
        clust.priorMean[FEH] = settings.cluster.priorMeans.Fe_H;
        clust.priorVar[FEH]  = settings.cluster.priorSigma.Fe_H;

        clust.mod = settings.cluster.starting.distMod;
        clust.priorMean[MOD] = settings.cluster.priorMeans.distMod;
        clust.priorVar[MOD]  = settings.cluster.priorSigma.distMod;

        clust.abs = settings.cluster.starting.Av;
        clust.priorMean[ABS] = fabs(settings.cluster.priorMeans.Av);
        clust.priorVar[ABS]  = settings.cluster.priorSigma.Av;

        clust.age = settings.cluster.starting.logAge;
        clust.priorMean[AGE] = settings.cluster.priorMeans.logAge;
        clust.priorVar[AGE]  = settings.cluster.priorSigma.logAge;

        clust.carbonicity = settings.cluster.starting.carbonicity;
        clust.priorMean[CARBONICITY] = settings.cluster.priorMeans.carbonicity;
        clust.priorVar[CARBONICITY]  = settings.cluster.priorSigma.carbonicity;

        clust.yyy = settings.cluster.starting.Y;
        clust.priorMean[YYY] = settings.cluster.priorMeans.Y;
        clust.priorVar[YYY]  = settings.cluster.priorSigma.Y;

        if (evoModels.IFMR <= 3)
        {
            clust.priorVar[IFMR_SLOPE] = 0.0;
            clust.priorVar[IFMR_INTERCEPT] = 0.0;
            clust.priorVar[IFMR_QUADCOEF] = 0.0;
        }
        else if (evoModels.IFMR <= 8)
        {
            clust.priorVar[IFMR_SLOPE] = 1.0;
            clust.priorVar[IFMR_INTERCEPT] = 1.0;
            clust.priorVar[IFMR_QUADCOEF] = 0.0;
        }
        else
        {
            clust.priorVar[IFMR_SLOPE] = 1.0;
            clust.priorVar[IFMR_INTERCEPT] = 1.0;
            clust.priorVar[IFMR_QUADCOEF] = 1.0;
        }

        /* set starting values for IFMR parameters */
        clust.ifmrSlope = clust.priorMean[IFMR_SLOPE] = 0.08;
        clust.ifmrIntercept = clust.priorMean[IFMR_INTERCEPT] = 0.65;

        if (evoModels.IFMR <= 10)
            clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.0001;
        else
            clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.08;

        clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

        for (auto &var : clust.priorVar)
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
    }

    void run();
    std::tuple<double, double, double> sampleMass(const Isochrone&, StellarSystem, int, int);
};


static double mirror(double var)
{
    // Mirror massRatio around 1 and 0
    if (var > 1.0)
    {
        var = 1.0 - (var - 1.0);

        return mirror(var);
    }
    else if (var < 0.0)
    {
        var = -var;

        return mirror(var);
    }

    return var;
}


// std::function<T(T)> propose
StellarSystem propose(std::mt19937 &gen, double massStepSize, double massRatioStepSize, double scale, StellarSystem s)
{
    double tMass = s.primary.mass + sampleT (gen, scale * massStepSize * massStepSize);
    double tMassRatio = s.getMassRatio() + sampleT (gen, scale * massRatioStepSize * massRatioStepSize);

    if (tMass < 0)
        tMass = -tMass;

    tMassRatio = mirror(tMassRatio);

    s.primary.mass = tMass;
    s.setMassRatio(tMassRatio);

    return s;
}


StellarSystem proposeCorrelated (std::mt19937 &gen, Matrix<double, 2, 2> const &propMatrix, StellarSystem s)
{
    array<double, 2> tDraws;

    for (auto &d : tDraws)
    {
        d = sampleT (gen, 1.0);
    }

    s.primary.mass += propMatrix.at(0).at(0) * tDraws[0];

    if (s.primary.mass < 0)
        s.primary.mass = -s.primary.mass;

    {
        double corrProps = 0;

        for (int k = 0; k < 2; ++k)
        {
            corrProps += propMatrix.at(1).at(k) * tDraws[k];
        }

        s.setMassRatio(mirror(s.getMassRatio() + corrProps));
    }

    return s;
}


// std::function<void(const T&)> checkPriors
void checkPriors(const Cluster &clust, const StellarSystem &s)
{
    if (s.getMassRatio() < 0.0)
        throw InvalidCluster("Low massRatio");
    else if (s.getMassRatio() > 1.0)
        throw InvalidCluster("High massRatio");
    else if (s.primary.mass <= 0.1)
        throw InvalidCluster("Low mass in primary");
    else if (s.primary.mass >= clust.getM_wd_up())
        throw InvalidCluster("High mass in primary");
}


std::tuple<double, double, double> Application::sampleMass(const Isochrone &isochrone, StellarSystem star, const int burnIters, const int iters)
{
    if (star.primary.mass < 0.1)
    {
        auto tRat = star.getMassRatio();

        star.primary.mass = 0.1;
        star.setMassRatio(tRat);
    }

    std::ofstream nullstream;

    // Proposal function bindings
    std::function<StellarSystem(StellarSystem)> burninProposal =
        std::bind(&propose, gen, massStepSize, massRatioStepSize, 10.0, _1);

    std::function<void(const StellarSystem&)> priorCheck = std::bind(&checkPriors, clust, _1);

    // This takes a special function binding because it's overridden.
    std::function<double(StellarSystem&)> logPost =
        std::bind<double (StellarSystem::*) (const Cluster&, const Model&, const Isochrone&) const>
        (&StellarSystem::logPost, _1, clust, evoModels, isochrone);

    // Start out with the burnin
    Chain<StellarSystem> burnin(static_cast<uint32_t>(std::uniform_int_distribution<>()(gen)), star, nullstream);
    burnin.run(burninProposal, logPost, priorCheck, burnIters);

    // Then do the main run
    if (settings.noBinaries)
    {
        Chain<StellarSystem> main(static_cast<uint32_t>(std::uniform_int_distribution<>()(gen)), burnin.get(), nullstream);
        main.run(burninProposal, logPost, priorCheck, iters);

        auto finalStar = main.get();

        try
        {
            double acceptedPosterior = exp(finalStar.logPost(clust, evoModels, isochrone));
            return std::make_tuple(finalStar.primary.mass, finalStar.getMassRatio(), acceptedPosterior);
        }
        catch (InvalidModelError &e)
        {
            double acceptedPosterior = -std::numeric_limits<double>::infinity();
            return std::make_tuple(99.999, 0.0, acceptedPosterior);
        }
    }
    else
    {
    // First, an abbreviated run to get a second covariance matrix
        try
        {
            std::function<StellarSystem(StellarSystem)> mainProposal =
                std::bind(&proposeCorrelated, gen, burnin.makeCholeskyDecomp(), _1);

            Chain<StellarSystem> main(static_cast<uint32_t>(std::uniform_int_distribution<>()(gen)), burnin.get(), nullstream);
            main.run(mainProposal, logPost, priorCheck, burnIters);

            // Then a secondary covariance matrix
            mainProposal = std::bind(&proposeCorrelated, gen, main.makeCholeskyDecomp(), _1);
            main.run(mainProposal, logPost, priorCheck, iters);

            auto finalStar = main.get();
            double acceptedPosterior = exp(finalStar.logPost(clust, evoModels, isochrone));

            return std::make_tuple(finalStar.primary.mass, finalStar.getMassRatio(), acceptedPosterior);
        }
        catch (NonPositiveDefiniteMatrix &e)
        {
            return std::make_tuple(99.99, 0.0, 0.0);
        }
    }
}


void Application::run()
{
    double fsLike;

    vector<StellarSystem> stars;

    {
        // Read photometry and calculate fsLike
        vector<double> filterPriorMin;
        vector<double> filterPriorMax;

        vector<string> filterNames;

        if (settings.run <= 0)
        {
            std::ifstream rData(settings.files.phot);

            if (!rData)
            {
                cerr << "***Error: Photometry file " << settings.files.phot << " was not found.***" << endl;
                cerr << ".at(Exiting...)" << endl;
                exit (1);
            }

            auto ret = base::utility::readPhotometry (rData, filterPriorMin, filterPriorMax, settings);

            filterNames = ret.first;
            stars = ret.second;
        }
        else
        {
            auto ret = base::utility::readPhotometryFromDB (filterPriorMin, filterPriorMax, settings);

            filterNames = ret.first;
            stars = ret.second;
        }

        evoModels.restrictFilters(filterNames);

        if (settings.cluster.index < 0 || static_cast<size_t>(settings.cluster.index) > filterNames.size())
        {
            cerr << "*** Error: " << settings.cluster.index << " is not a valid magnitude index.  Choose 0 through " << filterNames.size() - 1 << " ***"<< endl;
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

    auto sampledPars = readSampledParams (evoModels, settings);
    cout << "sampledPars[0].age = " << sampledPars.front().age << endl;

    /********** compile results *********/
    /*** now report sampled masses and parameters ***/

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

    for (size_t m = 0; m < sampledPars.size(); m++)
    {
        vector<std::pair<double, double>> masses;
        vector<double> memberships;

        clust.age = sampledPars.at(m).age;
        clust.feh = sampledPars.at(m).feh;
        clust.mod = sampledPars.at(m).distMod;
        clust.abs = sampledPars.at(m).abs;
        clust.yyy = sampledPars.at(m).y;
        clust.carbonicity = sampledPars.at(m).carbonicity;

        if (evoModels.IFMR >= 4)
        {
            clust.ifmrIntercept = sampledPars.at(m).ifmrIntercept;
            clust.ifmrSlope = sampledPars.at(m).ifmrSlope;
        }

        if (evoModels.IFMR >= 9)
        {
            clust.ifmrQuadCoef = sampledPars.at(m).ifmrQuadCoef;
        }

        /************ sample WD masses for different parameters ************/
        unique_ptr<Isochrone> isochrone(evoModels.mainSequenceEvol->deriveIsochrone(clust.feh, clust.yyy, clust.age));

        for (auto star : stars)
        {
            if (settings.noBinaries)
            {
                star.setMassRatio(0.0);
            }

            auto sampleTuple = sampleMass(*isochrone, star, settings.sampleMass.burnIters, settings.sampleMass.iters);

            double postClusterStar = std::get<2>(sampleTuple);
            postClusterStar *= (clust.getM_wd_up() - 0.15);

            masses.emplace_back(std::get<0>(sampleTuple), std::get<1>(sampleTuple));

            memberships.emplace_back(star.clustStarPriorDens * postClusterStar / (star.clustStarPriorDens * postClusterStar + (1.0 - star.clustStarPriorDens) * fsLike));
        }

        // massSampleFile << base::utility::format << sampledPars.at(m).age
        //                << base::utility::format << sampledPars.at(m).feh
        //                << base::utility::format << sampledPars.at(m).distMod
        //                << base::utility::format << sampledPars.at(m).abs;

        // if (evoModels.IFMR >= 4)
        // {
        //     massSampleFile << base::utility::format << sampledPars.at(m).ifmrIntercept
        //                    << base::utility::format << sampledPars.at(m).ifmrSlope;
        // }

        // if (evoModels.IFMR >= 9)
        // {
        //     massSampleFile << base::utility::format << sampledPars.at(m).ifmrQuadCoef;
        // }

        for (auto mass : masses)
        {
            massSampleFile << base::utility::format << mass.first
                           << base::utility::format << mass.second;
        }

        for (auto membership : memberships)
        {
            membershipFile << base::utility::format << membership;
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
    try
    {
        Settings settings;

        settings.loadSettings (argc, argv);

        if (settings.seed == std::numeric_limits<uint32_t>::max())
        {
            srand(std::time(0));
            settings.seed = rand();

            cout << "Seed: " << settings.seed << endl;
        }

        Application(settings).run();

        return 0;
    }
    catch (std::exception &e)
    {
        cerr << "\nException: " << e.what() << endl;
        return -1;
    }
}
