#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/format.hpp>

#include "Cluster.hpp"
#include "constants.hpp"
#include "Isochrone.hpp"
#include "Model.hpp"
#include "Settings.hpp"
#include "Star.hpp"

using std::cerr;
using std::endl;
using std::ofstream;
using std::string;
using std::vector;

class Application
{
  public:
    Application(Settings s)
        :  evoModels(makeModel(s)), settings(s)
    {}

    ~Application()
    {}

    void run();
    Isochrone interpolateIsochrone(const Cluster&);

  private:
    Model evoModels;
    const Settings settings;
};


void Application::run()
{
    vector<string> filters = evoModels.mainSequenceEvol->getAvailableFilters();

    for (size_t f = 0; f < filters.size(); ++f)
    {
        try
        {
            Filters::absCoeffs.at(filters.at(f));
        }
        catch(std::out_of_range &e)
        {
            filters.erase(filters.begin() + f);

            if (f > 0)
                f -= 1;
        }
    }

    evoModels.mainSequenceEvol->restrictToFilters(filters);

    Cluster clust;

    // Standard cluster parameters
    clust.feh = clust.priorMean[FEH] = settings.cluster.Fe_H;
    clust.mod = clust.priorMean[MOD] = settings.cluster.distMod;
    clust.abs = clust.priorMean[ABS] = fabs(settings.cluster.Av);
    clust.age = clust.priorMean[AGE] = settings.cluster.logClusAge;

    // IFMR parameters
    clust.ifmrSlope = clust.priorMean[IFMR_SLOPE] = 0.08;
    clust.ifmrIntercept = clust.priorMean[IFMR_INTERCEPT] = 0.65;


    // FilterSet-conditional parameters
    if (settings.whiteDwarf.wdModel == WdModel::MONTGOMERY)
    {
        clust.carbonicity = clust.priorMean[CARBONICITY] = settings.cluster.carbonicity;
    }
    else
    {
        clust.carbonicity = clust.priorMean[CARBONICITY] = 0.0;
    }

    if (settings.mainSequence.msRgbModel == MsModel::CHABHELIUM)
    {
        clust.yyy = clust.priorMean[YYY] = settings.cluster.Y;
    }
    else
    {
        clust.yyy = clust.priorMean[YYY] = 0.0;
    }

    if (evoModels.IFMR <= 10)
        clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.0001;
    else
        clust.ifmrQuadCoef = clust.priorMean[IFMR_QUADCOEF] = 0.08;

    clust.setM_wd_up(settings.whiteDwarf.M_wd_up);

    // Set Cluster's AGBT_zmass based on the model's agbTipMass
    clust.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(clust.feh, clust.yyy, clust.age);

    // Interpolate a new isochrone at 0.01M_0 steps
    auto interpIsochrone = interpolateIsochrone(clust);

    string filename = settings.files.output + ".cmd";

    ofstream fout(filename);

    if (!fout)
    {
        cerr << "***Error: Photometry file " << filename << " was not found.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    // Print file header
    {
        fout << boost::format("%10s") % "Mass";

        for (auto f : filters)
        {
            fout << boost::format(" %10s") % f;
        }

        fout << '\n';
    }

    for (auto eep : interpIsochrone.eeps)
    {
        fout << boost::format("%10.6f") % eep.mass;

        for (auto mag : eep.mags)
        {
            fout << boost::format(" %10.6f") % mag;
        }

        fout << '\n';
    }

    fout.close();
}


Isochrone Application::interpolateIsochrone(const Cluster &clust)
{
    Isochrone isochrone = evoModels.mainSequenceEvol->getIsochrone();

    vector<EvolutionaryPoint> eeps;;

    double mass = isochrone.eeps.front().mass;
    double delta = isochrone.agbTipMass() - mass;
    double deltaSteps = delta / 0.01;

    for ( int steps = 0; steps < deltaSteps; ++steps )
    {
        Star star;
        star.mass = mass + steps * 0.01;

        eeps.emplace_back(0, star.mass, star.getMags(clust, evoModels)); // 0 eep. Doesn't matter here.
    }

    return {clust.age, eeps};

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
