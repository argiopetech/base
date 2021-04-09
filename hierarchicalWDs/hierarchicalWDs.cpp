#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <getopt.h>
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
    unsigned int thin = 400;

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
        {"seed",      required_argument, 0,    0 },
        {"minLogAge", required_argument, 0,    1 },
        {"maxLogAge", required_argument, 0,    2 },
        {"samples",   required_argument, 0,    3 },
        {"thin",      required_argument, 0,    4 },
        {"help",      no_argument,       0, 0xFF },
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
                istringstream (string (optarg)) >> settings.thin;
                break;

            case 0xFF:
            case '?':
                cerr << "Usage: " << argv[0] << " [OPTION] FILES\n\n";
                cerr << "Fit a hierarchical age model to multiple white dwarfs given result files from singlePopMcmc. Outputs a single flat-text file with sampled gamma and tau squared (from model by Shijing Si et. al 2017) as first columns and one following column for each sampled age of each star in order of appearance on the command line.\n\n";
                cerr << "\t--help\n";
                cerr << "\t\tdisplay this help and exit\n\n";

                cerr << "\t--samples\n";
                cerr << "\t\tnumber of samples to take\n\n";
                cerr << "\t--thin\n";
                cerr << "\t\tnumber of iterations to run between samples\n\n";

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

    cout << "\nRunning simulation for " << settings.samples 
         << " samples thinning by " << settings.thin 
         << " (" << settings.samples * settings.thin << " total iterations.)\n";
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

vector<Result> sampleHierarchicalModel_FullyBayesian(vector<vector<double>> allAges, Settings settings)
{
    vector<Result> results;

    return results;
}

int main (int argc, char *argv[])
{
    Settings settings = loadCLISettings(argc, argv);

    reportSettings(settings);
    
    auto allAges = readResults(settings.resultFiles);

    auto results = sampleHierarchicalModel_FullyBayesian(allAges, settings);

    return 0;
}
