#ifndef WDEVOL_H
#define WDEVOL_H

#include <vector>

#include "Model.hpp"

double wdEvol (const Cluster &pCluster, const Model &, const std::vector<int>&, Star &pStar, int cmpnt);

#endif
