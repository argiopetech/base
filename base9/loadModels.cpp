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

void loadModels (Cluster *theCluster, const Model &evoModels, Settings &settings)
{
    setFilterNames (evoModels.filterSet);

    theCluster->carbonicity = settings.whiteDwarf.carbonicity;

    cout << "\nReading models..." << endl;

    if (evoModels.brownDwarfEvol == BARAFFE)
        loadBaraffe (settings.files.models);

    evoModels.mainSequenceEvol->loadModel(settings.files.models, settings.mainSequence.filterSet);

    loadWDCool (settings.files.models, evoModels.WDcooling);
    loadBergeron (settings.files.models, evoModels.filterSet);

    cout << "Models read." << endl;

}
