#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <vector>

#include "Chain.hpp"

#include "constants.hpp"
#include "densities.hpp"
#include "ifmr.hpp"
#include "Model.hpp"
#include "Settings.hpp"
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

struct ifmrGridControl
{
    array<double, NPARAMS> priorMean, priorVar;
    double minMag;
    double maxMag;
    int iStart;
    int modelSet;
    array<double, NPARAMS> start; /* starting points for grid evaluations */
    array<double, NPARAMS> end;   /* end points for grid evaluations */
};


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

/* Used in densities.c. */
array<double, NPARAMS> priorMean, priorVar;

/* TEMPORARY - global variable */
constexpr double const dMass1 = 0.0005;

Settings settings;


/*
 * read control parameters from input stream
 */
static void initIfmrGridControl (Chain *mc, Model &evoModels, struct ifmrGridControl *ctrl, Settings &s)
{
    ctrl->priorMean.at(CARBONICITY) = mc->clust.priorMean.at(CARBONICITY) = s.cluster.priorMeans.carbonicity;
    ctrl->priorVar.at(CARBONICITY) = s.cluster.priorSigma.carbonicity;

    ctrl->priorMean.at(FEH) = s.cluster.priorMeans.Fe_H;
    ctrl->priorVar.at(FEH)  = s.cluster.priorSigma.Fe_H;

    if (ctrl->priorVar.at(FEH) < 0.0)
    {
        ctrl->priorVar.at(FEH) = 0.0;
    }

    ctrl->priorMean.at(MOD) = s.cluster.priorMeans.distMod;
    ctrl->priorVar.at(MOD)  = s.cluster.priorSigma.distMod;
    if (ctrl->priorVar.at(MOD) < 0.0)
    {
        ctrl->priorVar.at(MOD) = 0.0;
    }

    ctrl->priorMean.at(ABS) = s.cluster.priorMeans.Av;
    ctrl->priorVar.at(ABS)  = s.cluster.priorSigma.Av;
    if (ctrl->priorVar.at(ABS) < 0.0)
    {
        ctrl->priorVar.at(ABS) = 0.0;
    }

    ctrl->priorMean.at(AGE) = s.cluster.priorMeans.logAge;
    ctrl->priorVar.at(AGE)  = s.cluster.priorSigma.logAge;

    ctrl->priorMean.at(YYY) = s.cluster.priorMeans.Y;
    ctrl->priorVar.at(YYY)  = s.cluster.priorSigma.Y;
    if (ctrl->priorVar.at(YYY) < 0.0)
    {
        ctrl->priorVar.at(YYY) = 0.0;
    }

    ctrl->priorMean.at(CARBONICITY) = s.cluster.priorMeans.carbonicity;
    ctrl->priorVar.at(CARBONICITY)  = s.cluster.priorSigma.carbonicity;
    if (ctrl->priorVar.at(CARBONICITY) < 0.0)
    {
        ctrl->priorVar.at(CARBONICITY) = 0.0;
    }

    ctrl->priorVar.at(IFMR_INTERCEPT) = 1.0;
    ctrl->priorVar.at(IFMR_SLOPE) = 1.0;
    if (evoModels.IFMR >= 9)
        ctrl->priorVar.at(IFMR_QUADCOEF) = 1.0;
    else
        ctrl->priorVar.at(IFMR_QUADCOEF) = 0.0;

    for (auto &var : ctrl->priorVar)
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

    // copy values to global variables
    priorVar.at(AGE) = ctrl->priorVar.at(AGE);
    priorVar.at(FEH) = ctrl->priorVar.at(FEH);
    priorVar.at(MOD) = ctrl->priorVar.at(MOD);
    priorVar.at(ABS) = ctrl->priorVar.at(ABS);
    priorVar.at(YYY) = ctrl->priorVar.at(YYY);
    priorVar.at(CARBONICITY) = ctrl->priorVar.at(CARBONICITY);
    priorVar.at(IFMR_INTERCEPT) = ctrl->priorVar.at(IFMR_INTERCEPT);
    priorVar.at(IFMR_SLOPE) = ctrl->priorVar.at(IFMR_SLOPE);
    priorVar.at(IFMR_QUADCOEF) = ctrl->priorVar.at(IFMR_QUADCOEF);

    priorMean.at(FEH) = ctrl->priorMean.at(FEH);
    priorMean.at(MOD) = ctrl->priorMean.at(MOD);
    priorMean.at(ABS) = ctrl->priorMean.at(ABS);
    priorMean.at(YYY) = ctrl->priorMean.at(YYY);
    priorMean.at(CARBONICITY) = ctrl->priorMean.at(CARBONICITY);

    /* prior values for linear IFMR */
    ctrl->priorMean.at(IFMR_SLOPE) = 0.08;
    ctrl->priorMean.at(IFMR_INTERCEPT) = 0.65;
    ctrl->priorMean.at(IFMR_QUADCOEF) = 0.0;
    priorMean.at(IFMR_SLOPE) = ctrl->priorMean.at(IFMR_SLOPE);
    priorMean.at(IFMR_INTERCEPT) = ctrl->priorMean.at(IFMR_INTERCEPT);
    priorMean.at(IFMR_QUADCOEF) = ctrl->priorMean.at(IFMR_QUADCOEF);

    /* open files for reading (data) and writing */

    ctrl->minMag = s.cluster.minMag;
    ctrl->maxMag = s.cluster.maxMag;

    ctrl->iStart = 0;

} // initIfmrGridControl


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
    mc->clust.feh = settings.cluster.starting.Fe_H;
    mc->clust.mod = settings.cluster.starting.Av;
    mc->clust.abs = settings.cluster.starting.Av;
    mc->clust.yyy = settings.cluster.starting.Y;
    mc->clust.age = settings.cluster.starting.logAge;
    mc->clust.carbonicity = settings.cluster.starting.carbonicity;
    mc->clust.ifmrIntercept = ctrl->priorMean.at(IFMR_INTERCEPT);
    mc->clust.ifmrSlope = ctrl->priorMean.at(IFMR_SLOPE);
    mc->clust.ifmrQuadCoef = ctrl->priorMean.at(IFMR_QUADCOEF);

    mc->clust.mean.at(AGE) = ctrl->priorMean.at(AGE);
    mc->clust.mean.at(YYY) = ctrl->priorMean.at(YYY);
    mc->clust.mean.at(MOD) = ctrl->priorMean.at(MOD);
    mc->clust.mean.at(FEH) = ctrl->priorMean.at(FEH);
    mc->clust.mean.at(ABS) = ctrl->priorMean.at(ABS);
    mc->clust.mean.at(CARBONICITY) = ctrl->priorMean.at(CARBONICITY);
    mc->clust.mean.at(IFMR_INTERCEPT) = ctrl->priorMean.at(IFMR_INTERCEPT);
    mc->clust.mean.at(IFMR_SLOPE) = ctrl->priorMean.at(IFMR_SLOPE);
    mc->clust.mean.at(IFMR_QUADCOEF) = ctrl->priorMean.at(IFMR_QUADCOEF);

    for (auto star : mc->stars)
    {
        star.clustStarProposalDens = star.clustStarPriorDens;   // Use prior prob of being clus star

        // find photometry for initial values of currentClust and mc->stars
        if (star.observedStatus == WD)
        {
            star.setMassRatio(0.0);
        }
    }
} // initChain

int main (int argc, char *argv[])
{
    Chain mc;
    struct ifmrGridControl ctrl;

    double fsLike;

    settings.loadSettings (argc, argv);

    if (settings.seed == std::numeric_limits<uint32_t>::max())
    {
        srand(std::time(0));
        settings.seed = rand();

        cout << "Seed: " << settings.seed << endl;
    }

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
            mc.stars = ret.second;
        }
        else
        {
            auto ret = base::utility::readPhotometryFromDB (filterPriorMin, filterPriorMax, settings);

            filterNames = ret.first;
            mc.stars = ret.second;
        }

        evoModels.restrictFilters(filterNames);

        if (   settings.cluster.index < 0
               || static_cast<size_t>(settings.cluster.index) > filterNames.size())
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

    mc.clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

    initChain (&mc, &ctrl);

    auto sampledPars = readSampledParams (evoModels, settings);
    cout << "sampledPars.at(0).age    = " << sampledPars.at(0).age << endl;
    cout << "sampledPars.at(last).age = " << sampledPars.back().age << endl;

    /* initialize WD logpost array and WD indices */
    double nWDLogPosts = (int) ceil ((mc.clust.getM_wd_up() - 0.15) / dMass1);

    /********** compile results *********/
    /*** now report sampled masses and parameters ***/

    // Open the file
    string filename = settings.files.output + ".wd";

    std::ofstream massFile;

    {
        string lf = filename + ".mass";

        massFile.open(lf);
        if (!massFile)
        {
            cerr << "***Error: File " << lf << " was not available for writing.***" << endl;
            cerr << ".at(Exiting...)" << endl;
            exit (1);
        }
    }

    std::ofstream membershipFile;

    {
        string lf = filename + ".membership";

        membershipFile.open(lf);
        if (!membershipFile)
        {
            cerr << "***Error: File " << lf << " was not available for writing.***" << endl;
            cerr << ".at(Exiting...)" << endl;
            exit (1);
        }
    }

    std::ofstream precLogAgeFile;

    if (settings.details)
    {
        string lf = filename + ".precLogAge";

        precLogAgeFile.open(lf);
        if (!precLogAgeFile)
        {
            cerr << "***Error: File " << lf << " was not available for writing.***" << endl;
            cerr << ".at(Exiting...)" << endl;
            exit (1);
        }
    }

    std::ofstream coolingAgeFile;

    if (settings.details)
    {
        string lf = filename + ".coolingAge";

        coolingAgeFile.open(lf);
        if (!coolingAgeFile)
        {
            cerr << "***Error: File " << lf << " was not available for writing.***" << endl;
            cerr << ".at(Exiting...)" << endl;
            exit (1);
        }
    }

    std::ofstream logTeffFile;

    if (settings.details)
    {
        string lf = filename + ".logTeff";

        logTeffFile.open(lf);
        if (!logTeffFile)
        {
            cerr << "***Error: File " << lf << " was not available for writing.***" << endl;
            cerr << ".at(Exiting...)" << endl;
            exit (1);
        }
    }

    std::ofstream loggFile;

    if (settings.details)
    {
        string lf = filename + ".logg";

        loggFile.open(lf);
        if (!loggFile)
        {
            cerr << "***Error: File " << lf << " was not available for writing.***" << endl;
            cerr << ".at(Exiting...)" << endl;
            exit (1);
        }
    }

    std::vector<double> us;
    us.resize(sampledPars.size() + 1);

    for (auto &u : us)
    {
        u = std::generate_canonical<double, 53>(gen);
    }

    for (size_t m = 0; m < sampledPars.size(); ++m)
    {
        Cluster internalCluster(mc.clust);
        std::vector<double> wdMass, clusMemPost;
        std::vector<double> wdLogPost;

        wdLogPost.resize(nWDLogPosts + 1);

        internalCluster.age = sampledPars.at(m).age;
        internalCluster.feh = sampledPars.at(m).feh;
        internalCluster.mod = sampledPars.at(m).distMod;
        internalCluster.abs = sampledPars.at(m).abs;

        if (evoModels.IFMR >= 4)
        {
            internalCluster.ifmrIntercept = sampledPars.at(m).ifmrIntercept;
            internalCluster.ifmrSlope = sampledPars.at(m).ifmrSlope;
        }

        if (evoModels.IFMR >= 9)
        {
            internalCluster.ifmrQuadCoef = sampledPars.at(m).ifmrQuadCoef;
        }

        /************ sample WD masses for different parameters ************/
        int iWD = 0;
        int im;
        double wdPostSum, maxWDLogPost, mass1;
        double postClusterStar;

        unique_ptr<Isochrone> isochrone(evoModels.mainSequenceEvol->deriveIsochrone(internalCluster.feh, internalCluster.yyy, internalCluster.age));

        for (auto star : mc.stars)
        {
            if (star.observedStatus == WD)
            {
                postClusterStar = 0.0;

                im = 0;

                for (mass1 = 0.15; mass1 < internalCluster.getM_wd_up(); mass1 += dMass1)
                {
                    /* condition on WD being cluster star */
                    star.primary.mass = mass1;
                    star.setMassRatio(0.0);

                    try
                    {
                        wdLogPost.at(im) = star.logPost (internalCluster, evoModels, *isochrone);
                        postClusterStar += exp (wdLogPost.at(im));
                    }
                    catch ( WDBoundsError &e )
                    {
                        cerr << e.what() << endl;

                        wdLogPost.at(im) = -HUGE_VAL;
                    }

                    im++;
                }
                im = 0;

                /* compute the maximum value */
                maxWDLogPost = wdLogPost.at(0);
                for (mass1 = 0.15; mass1 < internalCluster.getM_wd_up(); mass1 += dMass1)
                {
                    if (wdLogPost.at(im) > maxWDLogPost)
                        maxWDLogPost = wdLogPost.at(im);
                    im++;
                }

                /* compute the normalizing constant */
                wdPostSum = 0.0;
                im = 0;
                for (mass1 = 0.15; mass1 < internalCluster.getM_wd_up(); mass1 += dMass1)
                {
                    wdPostSum += exp (wdLogPost.at(im) - maxWDLogPost);
                    im++;
                }

                /* now sample a particular mass */
                double cumSum = 0.0;
                mass1 = 0.15;
                im = 0;
                while (cumSum < us.at(m) && mass1 < internalCluster.getM_wd_up())
                {
                    cumSum += exp (wdLogPost.at(im) - maxWDLogPost) / wdPostSum;
                    mass1 += dMass1;
                    im++;
                }
                mass1 -= dMass1;        /* maybe not necessary */

                wdMass.push_back(mass1);

                postClusterStar *= (internalCluster.getM_wd_up() - 0.15);

                clusMemPost.push_back(star.clustStarPriorDens * postClusterStar / (star.clustStarPriorDens * postClusterStar + (1.0 - star.clustStarPriorDens) * fsLike));
                iWD++;
            }
        }

        // massFile << base::utility::format << sampledPars.at(m).age
        //                << base::utility::format << sampledPars.at(m).feh
        //                << base::utility::format << sampledPars.at(m).distMod
        //                << base::utility::format << sampledPars.at(m).abs;

        // if (evoModels.IFMR >= 4)
        // {
        //     massFile << base::utility::format << sampledPars.at(m).ifmrIntercept
        //                    << base::utility::format << sampledPars.at(m).ifmrSlope;
        // }

        // if (evoModels.IFMR >= 9)
        // {
        //     massFile << base::utility::format << sampledPars.at(m).ifmrQuadCoef;
        // }

        for (size_t j = 0; j < wdMass.size(); j++)
        {
            auto precLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(internalCluster.feh, wdMass.at(j), internalCluster.yyy);

            double thisWDMass = intlFinalMassReln (internalCluster, evoModels, wdMass.at(j));
            double coolingAge = log10(exp10(internalCluster.age) - exp10(precLogAge));

            auto teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (internalCluster.age, internalCluster.carbonicity, precLogAge, thisWDMass);
            double logTeff = teffRadiusPair.first;

            auto logg = evoModels.WDAtmosphere->teffToLogg (logTeff, thisWDMass, WdAtmosphere::DA);

            massFile << base::utility::format << wdMass.at(j);
            membershipFile << base::utility::format << clusMemPost.at(j);

            if (settings.details)
            {
                precLogAgeFile << base::utility::format << precLogAge;
                coolingAgeFile << base::utility::format << coolingAge;
                logTeffFile << base::utility::format << logTeff;
                loggFile << base::utility::format << logg;
            }
        }

        massFile << endl;
        membershipFile << endl;

        if (settings.details)
        {
            precLogAgeFile << endl;
            coolingAgeFile << endl;
            logTeffFile << endl;
            loggFile << endl;
        };
    }

    massFile.close();
    membershipFile.close();

    if (settings.details)
    {
        precLogAgeFile.close();
        coolingAgeFile.close();
        logTeffFile.close();
        loggFile.close();
    }


    cout << "Part 2 completed successfully" << endl;

    return 0;
}
