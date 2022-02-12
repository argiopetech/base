#include <cstdlib>
#include <ctime>
#include <stdexcept>

#include "MpiMcmcApplication.hpp"
#include "IO/FieldStarLikelihood.hpp"
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


        SinglePopBackingStore *mcmcStore = nullptr;
        FieldStarLikelihoodBackingStore *fieldStarStore = nullptr;

        mcmcStore = new SinglePopMcmc_FileBackingStore(settings.files.output);
        fieldStarStore = new FieldStarLikelihood_FileBackingStore(settings.files.output);

        MpiMcmcApplication master(settings,
                                  std::move(mcmcStore),
                                  std::move(fieldStarStore));

        return master.run();
    }
    catch (exception &e)
    {
        cerr << "\nException: " << e.what() << endl;
        return -1;
    }
}
