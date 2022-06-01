#include <stdexcept>

#include "IO/FieldStarLikelihood.hpp"
#include "IO/MultiPopMcmc.hpp"
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

        mcmcStore      = new MultiPopMcmc_FileBackingStore(settings.files.output);
        paramsStore    = new StarParams_FileBackingStore(settings.files.output);

        // If you actually need FSLikelihood, uncomment the line below and comment the following line.
        // fieldStarStore = new FieldStarLikelihood_FileBackingStore(settings.files.output);
        fieldStarStore = new NullBackingStore<FieldStarLikelihoodRecord>();

        MpiMcmcApplication master(settings,
                                  std::move(mcmcStore),
                                  std::move(paramsStore),
                                  std::move(fieldStarStore));

        return master.run();
    }
    catch (exception &e)
    {
        cerr << "\nException: " << e.what() << endl;
        return -1;
    }

}
