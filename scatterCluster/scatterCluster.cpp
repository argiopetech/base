#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "Model.hpp"
#include "Settings.hpp"
#include "Star.hpp"
#include "Utility.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;
using std::setw;
using std::string;
using std::unordered_map;
using std::vector;

double getExposureTime(Settings &s, string f)
{
    try
    {
        return s.scatterCluster.exposures[f];
    }
    catch (std::out_of_range &e)
    {
        cerr << "No exposure field for filter \"" << f << "\" in configuration file" << endl;
        return 0.0;
    }
}

static unordered_map<string, std::pair<double, double>> s2nCoeffs_uncrowded = {
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


static unordered_map<string, std::pair<double, double>> s2nCoeffs_crowded = {
    {"U",       {9.33989,  0.3375778}},
    {"B",       {10.0478,  0.3462758}},
    {"V",       {10.48098, 0.3682010}},
    {"R",       {10.71151, 0.3837847}},
    {"I",       {10.61035, 0.3930941}},
    {"J",       {9.282385, 0.3862580}},
    {"H",       {9.197463, 0.3970419}},
    {"K",       {9.024068, 0.3985604}},
    {"UVf275w", {6.7302, 0.2090}},
    {"UVf336w", {6.9292, 0.2090}},
    {"UVf438w", {7.2452, 0.2090}},
    {"F606W",   {7.0797, 0.2090}},
    {"F814W",   {7.0797, 0.2090}}

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
double signalToNoise (const Settings &settings, double mag, double exptime, string filter)
{
    auto s2nCoeffs = settings.scatterCluster.crowded ? s2nCoeffs_crowded : s2nCoeffs_uncrowded;

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

bool meetsStageCutoff (const StellarSystem &s, const bool isFS)
{
    if (s.observedStatus == NSBH || s.observedStatus == DNE)
        return false;
    else if (s.observedStatus == WD && s.secondary.mass > 0.0 && !isFS)
        return false;                       // TEMPORARY KLUDGE -- ignore binaries of MS/RG + WDs and WD + WD
    else
        return true;
}


int main (int argc, char *argv[])
{
    Settings settings;

    vector<double> exposureTimes;

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

    if (settings.seed == std::numeric_limits<uint32_t>::max())
    {
        srand(std::time(0));
        settings.seed = rand();

        cout << "Seed: " << settings.seed << endl;
    }

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

        for (auto f : filters)
            exposureTimes.push_back(getExposureTime(settings, f));

        for (auto s : ret.second)
        {
            bool isFS = std::stoi(s.id) > 20000;
            if (meetsMagCutoff(settings, s) && meetsStageCutoff(s, isFS))
            {
                bool okay = true;

                for (size_t f = 0; f < filters.size(); ++f)
                {
                    auto filter = filters.at(f);
                    auto exposureTime = exposureTimes[f];

                    if (exposureTime > 0)
                    {
                        auto s2n    = signalToNoise (settings,
                                                     s.obsPhot.at(f),
                                                     exposureTime,
                                                     filter);

                        if (s2n < settings.scatterCluster.limitS2N)
                        {   // large photometric errors can lock mcmc during burnin
                            // Negative variances are ignored by mpiMcmc
                            // There's special code to check for negatives when we output stars (conversion to sigma)
                            okay = false;
                            break;
                        }
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
            auto s2n    = signalToNoise (settings,
                                         s.obsPhot.at(f),
                                         getExposureTime(settings, filter),
                                         filter);

            auto sigma  = 1 / s2n;

            auto targetSigma = settings.scatterCluster.crowded ? 0.04 : 0.01;

            if (sigma < targetSigma)
                sigma = targetSigma;

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

    // Output photometry header
    // All floating point have decimal precision = 3
    fout << setprecision(3) << fixed;

    fout << setw(4) << "id" << ' ';

    for (auto f : filters)
        if (getExposureTime(settings, f) > 0)
            fout << setw(6) << f << ' ';

    for (auto f : filters)
        if (getExposureTime(settings, f) > 0)
            fout << setw(6) << ("sig" + f) << ' ';

    fout << setw(7) <<     "mass1" << ' '
         << setw(9) << "massRatio" << ' '
         << setw(5) <<     "stage" << ' '
         << setw(7) <<   "Cmprior" << ' '
         << setw(6) <<    "useDBI" << ' '
         << endl;

    for (size_t i = 0; i < systems.size(); ++i)
    {
        auto s = systems.at(i);

        fout << setw(4) << (i + 1) << ' ';

        for (size_t f = 0; f < filters.size(); ++f)
        {
            if (exposureTimes[f] > 0)
            {
                fout << setw(6) << s.obsPhot[f]
                     << ' ';
            }
        }

        for (size_t f = 0; f < filters.size(); ++f)
        {
            if (exposureTimes[f] > 0)
            {
                auto v = s.variance[f];

                fout << setw(6) << (v < 0 ? v : sqrt (v))
                     << ' ';
            }
        }

        fout << setw(7) << s.primary.mass       << ' '
             << setw(9) << s.getMassRatio()     << ' '
             << setw(5) << s.observedStatus     << ' '
             << setw(7) << s.clustStarPriorDens << ' '
             << setw(6) << 1
             << endl;
    }


    return (0);
}
