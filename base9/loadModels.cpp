#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "evolve.hpp"
#include "loadModels.hpp"
#include "msRgbEvol.hpp"
#include "gBergMag.hpp"
#include "wdCooling.hpp"
#include "gBaraffeMag.hpp"
#include "FilterSet.hpp"

using std::cout;
using std::cerr;
using std::endl;

// Declared in parent program (mcmc, simCluster, makeCMD)
extern int verbose;

void loadModels (Cluster *theCluster, Model &evoModels, Settings &settings)
{
    setFilterNames (evoModels.filterSet);

    evoModels.IFMR = settings.whiteDwarf.ifmr;

    evoModels.WDcooling = settings.whiteDwarf.wdModel;

    if (evoModels.WDcooling < 0 || evoModels.WDcooling > 3)
    {
        cerr << "***Error: No models found for white dwarf filter set " << evoModels.WDcooling << ".***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    theCluster->carbonicity = settings.whiteDwarf.carbonicity;

    evoModels.brownDwarfEvol = settings.brownDwarf.bdModel;

    if (evoModels.brownDwarfEvol == BARAFFE)
        loadBaraffe (settings.files.models);

    evoModels.WDatm = BERGERON;

    cout << "\nReading models..." << endl;

    evoModels.mainSequenceEvol->loadModel(settings.files.models, settings.mainSequence.filterSet);

    loadWDCool (settings.files.models, evoModels.WDcooling);
    loadBergeron (settings.files.models, evoModels.filterSet);

    printf ("Models read.\n");

}
