#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
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
using std::ofstream;
using std::string;
using std::unique_ptr;
using std::vector;

class Application
{
  public:
    Application(Settings s)
        :  evoModels(makeModel(s)), wdEvoModels(makeModel(s)), settings(s)
    {}

    ~Application()
    {}

    void run();
    std::pair<Isochrone, Isochrone> interpolateIsochrone(const Cluster&, const Isochrone&);

  private:
    Model evoModels;
    Model wdEvoModels;

    const Settings settings;
};


void Application::run()
{
    Cluster clust;

    // Standard cluster parameters
    clust.feh = settings.cluster.starting.Fe_H;
    clust.priorMean[FEH] = settings.cluster.priorMeans.Fe_H;

    clust.mod = settings.cluster.starting.distMod;
    clust.priorMean[MOD] = settings.cluster.priorMeans.distMod;

    clust.abs = settings.cluster.starting.Av;
    clust.priorMean[ABS] = fabs(settings.cluster.priorMeans.Av);

    clust.age = settings.cluster.starting.logAge;
    clust.priorMean[AGE] = settings.cluster.priorMeans.logAge;

    clust.carbonicity = settings.cluster.starting.carbonicity;
    clust.priorMean[CARBONICITY] = settings.cluster.priorMeans.carbonicity;

    clust.yyy = settings.cluster.starting.Y;
    clust.priorMean[YYY] = settings.cluster.priorMeans.Y;

    // IFMR parameters
    clust.ifmrSlope = clust.priorMean[IFMR_SLOPE] = 0.08;
    clust.ifmrIntercept = clust.priorMean[IFMR_INTERCEPT] = 0.65;

    if (evoModels.IFMR <= 10)
        clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.0001;
    else
        clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.08;

    clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

    // Now setup the filters for MS and WDs
    {
        vector<string> msFilters = evoModels.mainSequenceEvol->getAvailableFilters();
        vector<string> filters;

        for (auto f : msFilters)
        {
            try
            {
                Filters::absCoeffs.at(f);

                filters.push_back(f);
            }
            catch(std::out_of_range &e)
            {}
        }

        evoModels.restrictFilters(filters);
    }

    {
        vector<string> msFilters = wdEvoModels.mainSequenceEvol->getAvailableFilters();
        vector<string> wdFilters = wdEvoModels.WDAtmosphere->getAvailableFilters();

        vector<string> filters;

        for ( auto m : msFilters )
        {
            for ( auto w : wdFilters )
            {
                if (m == w)
                {
                    try
                    {
                        // Throws if not available
                        Filters::absCoeffs.at(m);

                        // If we make it to this point, push
                        // This reminds me distinctly of goto...
                        filters.push_back(m);
                    }
                    catch(std::out_of_range &e)
                    {}

                    break;
                }
            }
        }

        wdEvoModels.restrictFilters(filters);
    }

    // Set Cluster's AGBT_zmass based on the model's agbTipMass
    unique_ptr<Isochrone> isochrone(evoModels.mainSequenceEvol->deriveIsochrone(clust.feh, clust.yyy, clust.age));

    // Interpolate a new isochrone at 0.01M_0 steps
    auto interpIsochrone = interpolateIsochrone(clust, *isochrone);

    {
        string filename = settings.files.output + ".ms";

        ofstream fout(filename);

        if (!fout)
        {
            cerr << "*** Error: Output file " << filename << " could not be opened.***" << endl;
            cerr << "*** [Exiting...]" << endl;
            exit (1);
        }

        // Print file header
        {
            fout << base::utility::format << "Mass";

            for (auto f : evoModels.mainSequenceEvol->getAvailableFilters())
            {
                fout << ' ' << base::utility::format << f;
            }

            fout << '\n';
        }

        for (auto eep : interpIsochrone.first.eeps)
        {
            fout << base::utility::format << eep.mass;

            for (auto mag : eep.mags)
            {
                fout << ' ' << base::utility::format << mag;
            }

            fout << '\n';
        }

        fout.close();
    }

    {
        string filename = settings.files.output + ".wd";

        ofstream fout(filename);

        if (!fout)
        {
            cerr << "*** Error: Output file " << filename << " could not be opened.***" << endl;
            cerr << "*** [Exiting...]" << endl;
            exit (1);
        }

        // Print file header
        {
            fout << base::utility::format << "Mass";

            for (auto f : wdEvoModels.WDAtmosphere->getAvailableFilters())
            {
                fout << ' ' << base::utility::format << f;
            }

            fout << '\n';
        }

        for (auto eep : interpIsochrone.second.eeps)
        {
            fout << base::utility::format << eep.mass;

            for (auto mag : eep.mags)
            {
                fout << ' ' << base::utility::format << mag;
            }

            fout << '\n';
        }

        fout.close();
    }
}


std::pair<Isochrone, Isochrone> Application::interpolateIsochrone(const Cluster &clust, const Isochrone &isochrone)
{
    vector<EvolutionaryPoint> eeps;
    vector<EvolutionaryPoint> wdEeps;

    {
        for ( size_t e = 0; e < isochrone.eeps.size() - 1; ++e)
        {
            double mass = isochrone.eeps.at(e).mass;
            double delta = isochrone.eeps.at(e + 1).mass - mass;
            double deltaSteps = delta / 80;

            for ( int steps = 0; steps < 80; ++steps )
            {
                StellarSystem star;
                star.primary.mass   = mass + deltaSteps * steps;
                star.secondary.mass = 0.0;

                eeps.emplace_back(e, star.primary.mass, star.deriveCombinedMags(clust, evoModels, isochrone));
            }
        }
    }

    {
        const double wdDelta = 0.01;
        double mass = isochrone.eeps.back().mass + wdDelta;

        while ( mass < clust.getM_wd_up() )
        {
            StellarSystem star;
            star.primary.mass   = mass;
            star.secondary.mass = 0.0;

            wdEeps.emplace_back(0, star.primary.mass, star.deriveCombinedMags(clust, wdEvoModels, isochrone));

            mass += wdDelta;
        }
    }

    return {{clust.age, eeps}, {clust.age, wdEeps}};

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
