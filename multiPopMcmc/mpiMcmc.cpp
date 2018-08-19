#include <stdexcept>

#include "IO/FieldStarLikelihood.hpp"
#include "IO/MultiPopMcmc.hpp"
#include "IO/Star.hpp"
#include "IO/StarParams.hpp"
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

        MultiPopBackingStore *mcmcStore = nullptr;
        StarParamsBackingStore *paramsStore = nullptr;
        FieldStarLikelihoodBackingStore *fieldStarStore = nullptr;
        StarBackingStore *photometryStore = nullptr;

        if (settings.files.backend == Backend::Sqlite)
        {
            auto tempStore = new MultiPopMcmc_SqlBackingStore(settings.files.output);
            mcmcStore      = tempStore;
            paramsStore    = new StarParams_SqlBackingStore(tempStore->runData());
            fieldStarStore = new FieldStarLikelihood_SqlBackingStore(tempStore->runData());
            photometryStore = new Star_SqlBackingStore(tempStore->runData());
        }
        else if (settings.files.backend == Backend::File)
        {
            mcmcStore      = new MultiPopMcmc_FileBackingStore(settings.files.output);
            paramsStore    = new StarParams_FileBackingStore(settings.files.output);
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
                                  std::move(paramsStore),
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
