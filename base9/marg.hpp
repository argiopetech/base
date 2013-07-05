#include <array>

#include "Model.hpp"


void setMags (double mag[][FILTS], int cmpnt, double *mass, Cluster &pCluster, Star &pStar, const Model&, std::array<double, 2>&);
double margEvolveWithBinary (Cluster &pCluster, Star &pStar, const Model&, std::array<double, 2>&);
