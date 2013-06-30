#include <string>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "msRgbEvol.hpp"

using std::string;

// void loadMSRgbModels (Cluster *pCluster, string path, int needFS)
// {

//     // The Yale models have issues creating field stars in simCluster, so if you need field stars
//     // and you are using the Yale models, also load the DSED models to create the field stars
//     if (pCluster->evoModels.mainSequenceEvol == YALE && needFS)
//         Dsed().loadMSRgbModels(path, pCluster->evoModels.filterSet);
// }
