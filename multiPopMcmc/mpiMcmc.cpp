#include <stdexcept>

#include "MpiMcmcApplication.hpp"

using std::cerr;
using std::endl;
using std::exception;

int main (int argc, char *argv[])
{
    cout << "Bayesian Analysis of Stellar Evolution (Multiple Populations)" << endl;

    Settings settings;

    // Setup settings
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

    MpiMcmcApplication master(settings);

    try
    {
        return master.run();
    }
    catch (exception &e)
    {
        cerr << "\nException: " << e.what() << endl;
        return -1;
    }
}
