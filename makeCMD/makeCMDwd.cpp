/* Author : Elliot Robinson
 *
 * Output photometry for expected WD stars between 1Mo and 7.5Mo
 */

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
        : wdEvoModels(makeModel(s)), settings(s)
    {}

    ~Application()
    {}

    void run();
    std::pair<Isochrone, Isochrone> interpolateIsochrone(const Cluster&, const Isochrone&);

  private:
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

    if (wdEvoModels.IFMR <= 10)
        clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.0001;
    else
        clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.08;

    clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

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

        wdEvoModels.restrictFilters(filters, false);
    }

    // Set Cluster's AGBT_zmass based on the model's agbTipMass
    unique_ptr<Isochrone> isochrone(wdEvoModels.mainSequenceEvol->deriveIsochrone(clust.feh, clust.yyy, clust.age));

    // Interpolate a new isochrone at 0.01M_0 steps
    auto interpIsochrone = interpolateIsochrone(clust, *isochrone);

    {
        string filename = settings.files.output + ".wda";

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
        string filename = settings.files.output + ".wdb";

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
    vector<EvolutionaryPoint> wDBEeps;
    vector<EvolutionaryPoint> wDAEeps;

    const double points[] = { 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,
                              4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5 };

    for ( auto p : points )
    {
        StellarSystem star;
        star.primary.mass   = p;
        star.secondary.mass = 0.0;

        star.primary.wdType = WdAtmosphere::DA;
        star.secondary.wdType = WdAtmosphere::DA;

        wDAEeps.emplace_back(0, star.primary.mass, star.deriveCombinedMags(clust, wdEvoModels, isochrone, settings.modIsParallax));

        star.primary.wdType = WdAtmosphere::DB;
        star.secondary.wdType = WdAtmosphere::DB;

        wDBEeps.emplace_back(0, star.primary.mass, star.deriveCombinedMags(clust, wdEvoModels, isochrone, settings.modIsParallax));
    }

    return {{clust.age, wDAEeps}, {clust.age, wDBEeps}};

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
