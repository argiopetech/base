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

using std::cerr;
using std::endl;

// Declared in parent program (mcmc, simCluster, makeCMD)
extern int verbose;

void loadModels (Cluster *theCluster, Model &evoModels, Settings &settings)
{
    char path[100] = "models/\0";

    evoModels.filterSet = settings.mainSequence.filterSet;

    if (evoModels.filterSet < 0 || evoModels.filterSet > 2)
    {
        printf ("***Error: No models found for filter set %d.***\n", evoModels.filterSet);
        printf ("[Exiting...]\n");
        exit (1);
    }

    setFilterNames (evoModels.filterSet);

    evoModels.IFMR = settings.whiteDwarf.ifmr;

    evoModels.WDcooling = settings.whiteDwarf.wdModel;

    if (evoModels.WDcooling < 0 || evoModels.WDcooling > 3)
    {
        printf ("***Error: No models found for white dwarf filter set %d.***\n", evoModels.WDcooling);
        printf ("[Exiting...]\n");
        exit (1);
    }

    theCluster->carbonicity = settings.whiteDwarf.carbonicity;

    evoModels.brownDwarfEvol = settings.brownDwarf.bdModel;

    if (evoModels.brownDwarfEvol == BARAFFE)
        loadBaraffe (path);

    evoModels.WDatm = BERGERON;

    printf ("\nReading models...\n");

    evoModels.mainSequenceEvol->loadModel(path, settings.mainSequence.filterSet);

    loadWDCool (path, evoModels.WDcooling);
    loadBergeron (path, evoModels.filterSet);

    printf ("Models read.\n");

}
