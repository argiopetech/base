#include <cstdlib>
#include <ctime>
#include <stdexcept>

#include "MpiMcmcApplication.hpp"
#include "IO/FieldStarLikelihood.hpp"
#include "IO/SinglePopMcmc.hpp"
#include "IO/Star.hpp"
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


        SinglePopBackingStore *mcmcStore = nullptr;
        FieldStarLikelihoodBackingStore *fieldStarStore = nullptr;
        StarBackingStore *photometryStore = nullptr;

        if (settings.files.backend == Backend::Sqlite)
        {
            auto tempStore = new SinglePopMcmc_SqlBackingStore(settings.files.output);
            mcmcStore = tempStore;
            fieldStarStore = new FieldStarLikelihood_SqlBackingStore(tempStore->runData());
            photometryStore = new Star_SqlBackingStore(tempStore->runData());
        }
        else if (settings.files.backend == Backend::File)
        {
            mcmcStore = new SinglePopMcmc_FileBackingStore(settings.files.output);
            fieldStarStore = new FieldStarLikelihood_FileBackingStore(settings.files.output);

            // We don't have a reason to save photometry that we've read from a file
            photometryStore = nullptr;
        }
        else
        {
            throw std::runtime_error("Invalid back end specified");
        }

        MpiMcmcApplication master(settings,
                                  std::move(mcmcStore),
                                  std::move(fieldStarStore),
                                  std::move(photometryStore));

        return master.run();
    }
    catch (exception &e)
    {
        cerr << "\nException: " << e.what() << endl;
        return -1;
    }
}
