#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/format.hpp>

#include "Model.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "Utility.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::unordered_map;
using std::vector;

double getExposureTime(Settings &s, string f)
{
    if (s.scatterCluster.exposures[f])
    {
        return s.scatterCluster.exposures[f];
    }
    else
    {
        cerr << "No exposure field for filter \"" << f << "\" in configuration file" << endl;
        return 0.0;
    }
}

static unordered_map<string, std::pair<double, double>> s2nCoeffs = {
    {"U",       {9.33989,  0.3375778}},
    {"B",       {10.0478,  0.3462758}},
    {"V",       {10.48098, 0.3682010}},
    {"R",       {10.71151, 0.3837847}},
    {"I",       {10.61035, 0.3930941}},
    {"J",       {9.282385, 0.3862580}},
    {"H",       {9.197463, 0.3970419}},
    {"K",       {9.024068, 0.3985604}},
    {"UVf275w", {7.804164, 0.2833815}},
    {"UVf336w", {7.779995, 0.2705339}},
    {"UVf438w", {7.958377, 0.2727768}},
    {"F606W",   {8.271906, 0.2647206}},
    {"F814W",   {8.272662, 0.2634004}}

    // {9.024068, 0.3985604},      // IRAC Blue
    // {9.024068, 0.3985604},      // IRAC Red
    // {9.024068, 0.3985604},      // Band 1
    // {9.024068, 0.3985604},      // Band 2
    // {9.024068, 0.3985604},      // Band 3
    // {9.024068, 0.3985604}       // Band 4
};

// This is an approximation to the results one would obtain in one
// hour with the KPNO 4m + Mosaic (UBVRI) or Flamingos (JHK) per band,
// assuming dark time, seeing=1.1", airmass=1.2, and then scaling from
// their by sqrt(exptime).  I further approximated the CCDTIME results
// (run on the NOAO webste) with linear fits of mag vs. log(S/N).
double signalToNoise (double mag, double exptime, string filter)
{
    double logS2N = s2nCoeffs[filter].first - s2nCoeffs[filter].second * mag;
    double s2n    = exp10 (logS2N);

    s2n *= sqrt (exptime);

    return s2n;
}

bool meetsMagCutoff (const Settings &settings, const StellarSystem &s)
{
    double mag;

    try
    {
        mag = s.obsPhot.at(settings.scatterCluster.relevantFilt);
    }
    catch (std::out_of_range &e)
    {
        cerr << "scatterCluster.relevantFilt out of range" << endl;
        throw;
    }

    if (s.observedStatus == BD)
        return true;
    else if (mag > 99.)
        return false;                       // check if real object.  if not, skip
    else if (mag < settings.scatterCluster.brightLimit)
        return false;                       // typically used to remove red giants
    else if (mag > settings.scatterCluster.faintLimit)
        return false;                       // typically used to remove WDs and lower MS
    else
        return true;
}

bool meetsStageCutoff (const StellarSystem &s)
{
    if (s.observedStatus == NSBH || s.observedStatus == NSBH)
        return false;
    else if (s.observedStatus == WD && s.secondary.mass > 0.0)
        return false;                       // TEMPORARY KLUDGE -- ignore binaries of MS/RG + WDs and WD + WD
    else
        return true;
}


int main (int argc, char *argv[])
{
    Settings settings;

    vector<string> filters;
    vector<StellarSystem> systems;

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

    {
        // Read photometry
        vector<double> filterPriorMin;
        vector<double> filterPriorMax;

        // open files for reading (data) and writing
        // rData implcitly relies on going out of scope to close the photometry file
        // This is awful, but pretty (since this code is, at time of writing, in restricted, anonymous scope)
        string filename = settings.files.output + ".sim.out";
        std::ifstream photin(filename);

        if (!photin)
        {
            cerr << "***Error: Photometry file " << filename << " was not found.***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (-1);
        }

        auto ret = base::utility::readPhotometry (photin, filterPriorMin, filterPriorMax, settings);

        filters = ret.first;

        for (auto s : ret.second)
        {
            if (meetsMagCutoff(settings, s) && meetsStageCutoff(s))
            {
                bool okay = true;

                for (size_t f = 0; f < filters.size(); ++f)
                {
                    auto filter = filters.at(f);
                    auto s2n    = signalToNoise (s.obsPhot.at(f),
                                                 getExposureTime(settings, filter),
                                                 filter);

                    if (s2n < settings.scatterCluster.limitS2N)
                    {   // large photometric errors can lock mcmc during burnin
                        // Negative variances are ignored by mpiMcmc
                        // There's special code to check for negatives when we output stars (conversion to sigma)
                        okay = false;
                        break;
                    }
                }

                if (okay)
                {
                    systems.push_back(s);
                }
            }
        }

        if (systems.empty())
        {
            cerr << "No valid input photometry after considering mag and stage cuttoffs" << endl;
            exit(0);
        }
    }

    std::mt19937 gen(uint32_t(settings.seed * uint32_t(2654435761))); // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)

    for (auto &s : systems)
    {
        for (size_t f = 0; f < filters.size(); ++f)
        {
            auto filter = filters.at(f);
            auto s2n    = signalToNoise (s.obsPhot.at(f),
                                         getExposureTime(settings, filter),
                                         filter);

            auto sigma  = 1 / s2n;

            if (sigma < 0.01)
                sigma = 0.01;

            s.variance.at(f) = sigma * sigma;

            std::normal_distribution<double> normDist(0., sigma);
            s.obsPhot.at(f) += normDist(gen);
        }
    }

    string filename = settings.files.output + ".scatter.out";
    std::ofstream fout(filename);

    if (!fout)
    {
        cerr << "***Error: Photometry file " << filename << " was not found.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (-1);
    }

    //Output photometry header
    fout << boost::format("%4s ") % "id";

    for (auto f : filters)
        fout << boost::format("%6s ") % f;

    for (auto f : filters)
        fout << boost::format("%6s ") % ("sig" + f);

    fout << boost::format("%7s %9s %5s %7s %6s\n") % "mass1" % "massRatio" % "stage" % "Cmprior" % "useDBI";

    for (size_t i = 0; i < systems.size(); ++i)
    {
        auto s = systems.at(i);

        fout << boost::format("%4d ") % (i + 1);

        for (auto f : s.obsPhot)
            fout << boost::format("%6.3f ") % f;

        for (auto f : s.variance)
        {
            fout << boost::format("%6.3f ") % (f < 0
                                                 ? f
                                                 : sqrt (f));
        }

        fout << boost::format("%7.3f %9.3f %5d %7.3f %6d\n")
            % s.primary.mass
            % s.getMassRatio()
            % s.observedStatus
            % s.clustStarPriorDens
            % 1;
    }


    return (0);
}
