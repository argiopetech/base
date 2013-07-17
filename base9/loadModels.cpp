#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "evolve.hpp"
#include "loadModels.hpp"
#include "gBergMag.hpp"
#include "WdCoolingModel.hpp"
#include "gBaraffeMag.hpp"

using std::cout;
using std::cerr;
using std::endl;

void loadModels (Cluster &theCluster, const Model &evoModels, const Settings &settings)
{
    theCluster.carbonicity = settings.whiteDwarf.carbonicity;

}
