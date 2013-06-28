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

// Declared in parent program (mcmc, simCluster, makeCMD)
extern int verbose;

void loadModels (int needFS, Cluster *theCluster, Settings &settings)
{
    char path[100] = "models/\0";

    theCluster->evoModels.mainSequenceEvol = settings.mainSequence.msRgbModel;

    if (theCluster->evoModels.mainSequenceEvol < 0 || theCluster->evoModels.mainSequenceEvol > 3)
    {
        printf ("***Error: No models found for main sequence evolution model %d.***\n", theCluster->evoModels.mainSequenceEvol);
        printf ("[Exiting...]\n");
        exit (1);
    }

    theCluster->evoModels.filterSet = settings.mainSequence.filterSet;

    if (theCluster->evoModels.filterSet < 0 || theCluster->evoModels.filterSet > 2)
    {
        printf ("***Error: No models found for filter set %d.***\n", theCluster->evoModels.filterSet);
        printf ("[Exiting...]\n");
        exit (1);
    }

    setFilterNames (theCluster->evoModels.filterSet);

    theCluster->evoModels.IFMR = settings.whiteDwarf.ifmr;

    theCluster->evoModels.WDcooling = settings.whiteDwarf.wdModel;

    if (theCluster->evoModels.WDcooling < 0 || theCluster->evoModels.WDcooling > 3)
    {
        printf ("***Error: No models found for white dwarf filter set %d.***\n", theCluster->evoModels.WDcooling);
        printf ("[Exiting...]\n");
        exit (1);
    }

    theCluster->carbonicity = settings.whiteDwarf.carbonicity;

    theCluster->evoModels.brownDwarfEvol = settings.brownDwarf.bdModel;

    if (theCluster->evoModels.brownDwarfEvol == BARAFFE)
        loadBaraffe (path);

    theCluster->evoModels.WDatm = BERGERON;

    printf ("\nReading models...\n");

    loadMSRgbModels (theCluster, path, needFS);
    loadWDCool (path, theCluster->evoModels.WDcooling);
    loadBergeron (path, theCluster->evoModels.filterSet);

    printf ("Models read.\n");

}
