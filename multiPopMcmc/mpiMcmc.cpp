#include <stdexcept>

#include "IO/MultiPopMcmc.cpp"
#include "IO/StarParams.cpp"
#include "MpiMcmcApplication.hpp"

using std::cerr;
using std::endl;
using std::exception;


int main (int argc, char *argv[])
{
    Settings settings;

    try
    {
        // Setup settings
        settings.loadSettings (argc, argv);

        cout << "Bayesian Analysis of Stellar Evolution (Multiple Populations)" << endl;

        if (settings.seed == std::numeric_limits<uint32_t>::max())
        {
            srand(std::time(0));
            settings.seed = rand();

            cout << "Seed: " << settings.seed << endl;
        }

        // TODO - Make this read settings
        auto mcmcStore   = new MultiPopMcmc_FileBackingStore(settings.files.output + ".res");
        auto paramsStore = new StarParams_FileBackingStore(settings.files.output + ".starParams");

        MpiMcmcApplication master(settings, std::move(mcmcStore), std::move(paramsStore));

        return master.run();
    }
    catch (exception &e)
    {
        cerr << "\nException: " << e.what() << endl;
        return -1;
    }

}
