/* Author: Elliot Robinson
 *
 * Outputs CMDs for every line in a .res output file
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

#include "Cluster.hpp"
#include "constants.hpp"
#include "Isochrone.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "Utility.hpp"

using std::cerr;
using std::endl;
using std::fabs;
using std::stringstream;
using std::ofstream;
using std::string;
using std::unique_ptr;
using std::vector;

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

static vector<clustPar> readSampledParams (Model &evoModels, const Settings &s)
{
    string line;

    vector<clustPar> sampledPars;

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

    return sampledPars;
}


class Application
{
  public:
    Application(Settings s)
        :  evoModels(makeModel(s)), settings(s)
    {}

    ~Application()
    {}

    void run();
    Isochrone interpolateIsochrone(const Cluster&, const Isochrone&);
    Isochrone interpolateWDIsochrone(const Cluster&, const Isochrone&);

  private:
    Model evoModels;
    const Settings settings;
};


void Application::run()
{

    // for (size_t f = 0; f < filters.size(); ++f)
    // {
    //     try
    //     {
    //         Filters::absCoeffs.at(filters.at(f));
    //     }
    //     catch(std::out_of_range &e)
    //     {
    //         filters.erase(filters.begin() + f);

    //         if (f > 0)
    //             f -= 1;
    //     }
    // }

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

        evoModels.restrictFilters(filterNames);
    }

    // evoModels.restrictFilters(filters);
    vector<string> filters = evoModels.mainSequenceEvol->getAvailableFilters();

    Cluster clust;

    auto sampledParams = readSampledParams (evoModels, settings);

    double sampledMin = std::numeric_limits<double>::max();
    double sampledMax = std::numeric_limits<double>::lowest();

    for (auto s : sampledParams)
    {
        sampledMin = std::min(sampledMin, s.logPost);
        sampledMax = std::max(sampledMax, s.logPost);
    }

    auto alpha = [=](double x){
        double delta = sampledMax - sampledMin;
        return (x - sampledMin) / delta;
    };


    string filename = settings.files.output + ".sampleCmds";

    ofstream fout(filename);

    if (!fout)
    {
        cerr << "*** Error: Output file " << filename << " could not be opened.***" << endl;
        cerr << "*** [Exiting...]" << endl;
        exit (1);
    }

    // Print file header
    {
        for (auto f : filters)
        {
            fout << base::utility::format << f << ' ';
        }

        fout << base::utility::format << "alpha";

        fout << '\n';
    }

    for (auto s : sampledParams)
    {
        // Standard cluster parameters
        clust.feh = clust.priorMean[FEH] = s.feh;
        clust.mod = clust.priorMean[MOD] = s.distMod;
        clust.abs = clust.priorMean[ABS] = fabs(s.abs);
        clust.age = clust.priorMean[AGE] = s.age;

        // IFMR parameters
        clust.ifmrSlope = clust.priorMean[IFMR_SLOPE] = 0.08;
        clust.ifmrIntercept = clust.priorMean[IFMR_INTERCEPT] = 0.65;

        // FilterSet-conditional parameters
        clust.carbonicity = clust.priorMean[CARBONICITY] = s.carbonicity;
        clust.yyy = clust.priorMean[YYY] = s.y;

        if (evoModels.IFMR <= 10)
            clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.0001;
        else
            clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.08;

        clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

        // Set Cluster's AGBT_zmass based on the model's agbTipMass
        unique_ptr<Isochrone> isochrone(evoModels.mainSequenceEvol->deriveIsochrone(clust.feh, clust.yyy, clust.age));

        // Interpolate a new isochrone
        auto interpIsochrone = interpolateIsochrone(clust, *isochrone);

        for (auto eep : interpIsochrone.eeps)
        {
            for (auto mag : eep.mags)
            {
                fout << base::utility::format << mag << ' ';
            }

            fout << base::utility::format << alpha(s.logPost);

            fout << '\n';
        }

    }

    fout.close();
}


Isochrone Application::interpolateIsochrone(const Cluster &clust, const Isochrone &isochrone)
{
    vector<EvolutionaryPoint> eeps;

    for ( size_t e = 0; e < isochrone.eeps.size(); ++e)
    {
        double mass = isochrone.eeps.at(e).mass;

        if (mass > isochrone.agbTipMass())
            continue;

        {
            StellarSystem star;
            star.primary.mass   = mass;
            star.secondary.mass = 0.0;

            eeps.emplace_back(e, star.primary.mass, star.deriveCombinedMags(clust, evoModels, isochrone));
        }
    }

    return {clust.age, eeps};
}


Isochrone Application::interpolateWDIsochrone(const Cluster &clust, const Isochrone &isochrone)
{
    vector<EvolutionaryPoint> eeps;

    double mass = 0;
    int e = 0;
    do
    {
        mass = isochrone.agbTipMass() + ++e * 0.05;

        if (mass > 8.0)
            continue;

        {
            StellarSystem star;
            star.primary.mass   = mass;
            star.secondary.mass = 0.0;

            eeps.emplace_back(e, star.primary.mass, star.deriveCombinedMags(clust, evoModels, isochrone));
        }
    } while (mass < 5.0);

    return {clust.age, eeps};
}


int main (int argc, char *argv[])
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
