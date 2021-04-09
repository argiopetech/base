#include <cstdio>
#include <cstdlib>
#include <ctime>
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
using std::vector;

struct Settings
{
    unsigned int seed = std::numeric_limits<uint32_t>::max();
    vector<string> resultFiles;
};

Settings loadCLISettings(int argc, char *argv[])
{
    int c;
    int option_index = 0;

    Settings settings;

    static struct option long_options[] = {
        {"seed",    required_argument, 0,  0 },
        {"help",    no_argument,       0,  1 },
        {0,         0,                 0,  0 }
    };

    while ((c = getopt_long(argc, argv, "s:", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
            case 's':
                istringstream (string (optarg)) >> settings.seed;
                break;

            case 1:
            case '?':
                cerr << "Usage: " << argv[0] << " [OPTION] FILES\n";
                cerr << "Run a \n\n";
                cerr << "\t--help\n";
                cerr << "\t\tdisplay this help and exit\n";
                cerr << "\t-s, --seed\n";
                cerr << "\t\tSet the seed for the random number generator\n";
                exit (EXIT_FAILURE);
                break;

            default:
                abort();
        }
    }

    if (optind - argc < 2) {
        while (optind < argc)
        {
            settings.resultFiles.push_back(argv[optind++]);
        }
    }
    else
    {
        cerr << "Must pass at least two result files. See --help for details.";
            
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
    cout << "Seed: " << settings.seed << endl;

    cout << "Loading results from:\n";

    for (auto f : settings.resultFiles)
    {
        cout << "\t- " << f << "\n";
    }
}

int main (int argc, char *argv[])
{
    Settings settings = loadCLISettings(argc, argv);

    reportSettings(settings);
    
    return 0;
}
