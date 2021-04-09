#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::istringstream;
using std::string;
using std::stringstream;
using std::vector;

struct Settings
{
    unsigned int seed = std::numeric_limits<uint32_t>::max();

    unsigned int samples = 10000;
    unsigned int chainLength = 400;

    double minLogAge = 9.5;
    double maxLogAge = 10.176;

    vector<string> resultFiles;
};

Settings loadCLISettings(int argc, char *argv[])
{
    int c;
    int option_index = 0;

    Settings settings;

    static struct option long_options[] = {
        {"seed",        required_argument, 0,    0 },
        {"minLogAge",   required_argument, 0,    1 },
        {"maxLogAge",   required_argument, 0,    2 },
        {"samples",     required_argument, 0,    3 },
        {"chainLength", required_argument, 0,    4 },
        {"help",        no_argument,       0, 0xFF },
        {0, 0, 0, 0 }
    };

    while ((c = getopt_long(argc, argv, "s:", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
            case 's':
                istringstream (string (optarg)) >> settings.seed;
                break;

            case 1:
                istringstream (string (optarg)) >> settings.minLogAge;
                break;

            case 2:
                istringstream (string (optarg)) >> settings.maxLogAge;
                break;

            case 3:
                istringstream (string (optarg)) >> settings.samples;
                break;

            case 4:
                istringstream (string (optarg)) >> settings.chainLength;
                break;

            case 0xFF:
            case '?':
                cerr << "Usage: " << argv[0] << " [OPTION] FILES\n\n";
                cerr << "Fit a hierarchical age model to multiple white dwarfs given result files from singlePopMcmc. Outputs a single flat-text file with sampled gamma and tau squared (from model by Shijing Si et. al 2017) as the first two columns and one following column for each sampled age of each star in order of appearance on the command line.\n\n";
                cerr << "\t--help\n";
                cerr << "\t\tdisplay this help and exit\n\n";

                cerr << "\t--samples\n";
                cerr << "\t\tnumber of samples to take\n\n";
                cerr << "\t--chainLength\n";
                cerr << "\t\tnumber of iterations to run for each sample\n\n";

                cerr << "\t--maxLogAge\n";
                cerr << "\t\tmaximum age for the gamma model parameter\n\n";
                cerr << "\t--minLogAge\n";
                cerr << "\t\tminimum age for the gamma model parameter\n\n";

                cerr << "\t-s, --seed\n";
                cerr << "\t\tSet the seed for the random number generator\n";

                exit (EXIT_FAILURE);
                break;

            default:
                abort();
        }
    }

    if ((argc - optind) >= 2) {
        while (optind < argc)
        {
            settings.resultFiles.push_back(argv[optind++]);
        }
    }
    else
    {
        cerr << "Must pass at least two result files. See --help for details.";

        exit (EXIT_FAILURE);
    }

    if (settings.seed == std::numeric_limits<uint32_t>::max())
    {
        srand(std::time(0));

        // Applies Knuth's multiplicative hash for obfuscation (TAOCP Vol. 3)
        settings.seed = rand() * uint32_t(2654435761);
    }

    return settings;
}

void reportSettings(Settings settings)
{
    cout << "Seed: " << settings.seed << "\n\n";

    cout << "Minimum age: " << settings.minLogAge << " log years\n";
    cout << "Maximum age: " << settings.maxLogAge << " log years\n\n";

    cout << "Loading results from:\n";

    for (auto f : settings.resultFiles)
    {
        cout << "\t- " << f << "\n";
    }

    cout << "\nTaking " << settings.samples
         << " samples with " << settings.chainLength
         << " iteration chains (" << settings.samples * settings.chainLength << " total iterations.)\n";
}

static vector<vector<double>> readResults (vector<string> resultFiles)
{
    vector<vector<double>> sampledAges;

    string line;

    for (string file : resultFiles)
    {
        sampledAges.emplace_back();

        std::ifstream parsFile;
        parsFile.open(file);

        getline(parsFile, line); // Skip header

        while (getline(parsFile, line))
        {
            stringstream in(line);

            double newAge;
            in >> newAge;

            // Only take stage 3 (main run) iterations
            if (line.back() == '3')
                sampledAges.back().push_back(newAge);
        }

        parsFile.close();
    }

    return sampledAges;
}

struct Result
{
    Result(double gamma, double tauSquared, vector<double> sampledAges)
        : gamma(gamma), tauSquared(tauSquared), sampledAges(sampledAges)
    {;}

    double gamma;
    double tauSquared;
    vector<double> sampledAges;
};

vector<double> sampleStarAges(vector<vector<double>> const &allAges, std::mt19937 &gen)
{
    vector<double> sampledAges;

    for (auto star = 0u; star < allAges.size(); ++star)
    {
        std::uniform_int_distribution<size_t> randomIndex(0, allAges[star].size() - 1);

        sampledAges.push_back(allAges[star][randomIndex(gen)]);
    }

    return sampledAges;
}

constexpr double sqr(double a)
{
    return a * a;
}

double meanOf(vector<double> const &values)
{
    double mean = 0;

    for (auto v : values)
    {
        mean += v;
    }

    return mean / values.size();
}

double totalSumOfSquaresOf(vector<double> const &values, double mean)
{
    double result = 0;

    for (auto v : values)
    {
        result += sqr(v - mean);
    }

    return result;
}

double truncatedNormallyDistributedDraw(double mean, double standardDeviation, double min, double max, std::mt19937 &gen)
{
    int count = 0;
    bool once = false;

    double result;

    std::normal_distribution<double> normalDistr(mean, standardDeviation);

    do
    {
        result = normalDistr(gen);
        count += 1;

        if (!once && count > 100000)
        {
            once = true;

            cerr << "\rStuck drawing a truncated, normally distributed value. Perhaps your min or max bound needs adjustment?\n";
        }
    } while ((result < min) || (result > max));

    return result;
}

double normalPDF(double mean, double standardDeviation, double point)
{
    double left = 1 / (standardDeviation * sqrt(2 * M_PI));
    double right = exp(-0.5 * sqr((point - mean) / standardDeviation));

    return left * right;
}

void percentComplete(int percent)
{
    cout << '\r' << percent << "% complete" << std::flush;
}

int reportPercentage(int lastPercent, int current, int target)
{
    int percent = (100 * current) / target;

    if (percent > lastPercent)
    {
        percentComplete(percent);

        lastPercent = percent;
    }

    return percent;
}

vector<Result> sampleHierarchicalModel_FullyBayesian(vector<vector<double>> const &allAges, Settings const settings, std::mt19937 &gen)
{
    vector<Result> results;

    int lastPercent = 0;
    percentComplete(0);

    auto nStars = allAges.size();

    // Starting values
    auto sampledAges = sampleStarAges(allAges, gen);

    double tauSquaredScale = totalSumOfSquaresOf(sampledAges, meanOf(sampledAges));
    double tauSquared      = tauSquaredScale / (nStars - 1);
    double gamma           = truncatedNormallyDistributedDraw(meanOf(sampledAges), sqrt(tauSquared / nStars), settings.minLogAge, settings.maxLogAge, gen);

    results.emplace_back(gamma, tauSquared, sampledAges);

    for (auto sample = 0u; sample < settings.samples; ++sample)
    {
        double tau = sqrt(tauSquared);

        for (auto i = 0u; i < settings.chainLength; ++i)
        {
            auto newAges = sampleStarAges(allAges, gen);

            for (auto star = 0u; star < nStars; ++star)
            {
                auto probabilityRatio = normalPDF(gamma, tau, newAges[star]) /
                                     // ---------------------------------------
                                        normalPDF(gamma, tau, sampledAges[star]);

                double u = std::generate_canonical<double, 53>(gen);

                if (probabilityRatio <= u)
                {
                    sampledAges[star] = newAges[star];
                }
            }
        }

        tauSquaredScale = totalSumOfSquaresOf(sampledAges, gamma);

        std::gamma_distribution<double> gammaDist((nStars - 1) / 2, tauSquaredScale / 2);

        tauSquared = 1 / gammaDist(gen); // Inverts the draw from the gamma distribution, giving inverse-gamma
        gamma      = truncatedNormallyDistributedDraw(meanOf(sampledAges), sqrt(tauSquared / nStars), settings.minLogAge, settings.maxLogAge, gen);

        results.emplace_back(gamma, tauSquared, sampledAges);

        lastPercent = reportPercentage(lastPercent, sample, settings.samples);
   }

    percentComplete(100);

    return results;
}

void exportResults(vector<Result> const &results)
{
    auto format = [](std::ostream& out) -> std::ostream&
        {
            return out << std::setw(11) << std::fixed << std::setprecision(6);
        };

    string filename = "hierarchicalWDs.out";

    std::ofstream fout(filename);

    if(!fout)
    {
        throw std::runtime_error(filename + " was not available for writing.");
    }

    for (auto result : results)
    {
        fout << format << result.gamma 
             << ' ' << format << result.tauSquared;

        for (auto age : result.sampledAges)
        {
            fout << ' ' << format << age;
        }

        fout << "\n";
    }

    fout.close();
}

int main (int argc, char *argv[])
{
    Settings settings = loadCLISettings(argc, argv);

    reportSettings(settings);

    auto allAges = readResults(settings.resultFiles);

    std::mt19937 gen(settings.seed);
    gen.discard(10000); // Warms up the PRNG

    auto results = sampleHierarchicalModel_FullyBayesian(allAges, settings, gen);

    exportResults(results);

    return 0;
}
