#include <cstdlib>
#include <ctime>
#include <stdexcept>

#include "MpiMcmcApplication.hpp"
#include "IO/SinglePopMcmc.hpp"
#include "IO/Records.hpp"

using std::cerr;
using std::endl;
using std::exception;


int main (int argc, char *argv[])
{
    try
    {
        cout << "Bayesian Analysis of Stellar Evolution" << endl;

        Settings settings;

        // Setup settings
        settings.loadSettings (argc, argv);

        if (settings.seed == std::numeric_limits<uint32_t>::max())
        {
            srand(std::time(0));
            settings.seed = rand();

            cout << "Seed: " << settings.seed << endl;
        }


        SinglePopBackingStore *store = nullptr;

        if (settings.files.backend == Backend::Sqlite)
        {
            store = new SinglePopMcmc_SqlBackingStore(settings.files.output);
        }
        else if (settings.files.backend == Backend::File)
        {
            store = new SinglePopMcmc_FileBackingStore(settings.files.output);
        }
        else
        {
            throw std::runtime_error("Invalid back end specified");
        }

        MpiMcmcApplication master(settings, std::move(store));

        return master.run();
    }
    catch (exception &e)
    {
        cerr << "\nException: " << e.what() << endl;
        return -1;
    }
}
