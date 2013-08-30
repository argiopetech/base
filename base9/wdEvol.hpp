#ifndef WDEVOL_H
#define WDEVOL_H

#include <array>
#include <vector>

#include "Model.hpp"

double wdEvol (const Cluster &pCluster, const Model &, const std::vector<int>&, std::array<double, FILTS>&, Star &pStar, int cmpnt);

#endif
