#include <array>

#include "Model.hpp"

void deriveCombinedMags (double mag[][FILTS], double clusterAv, double *flux, Cluster &pCluster, Star &pStar, const Model &evoModels, const std::vector<int>&);
void setMags (double mag[][FILTS], int cmpnt, double *mass, Cluster &pCluster, Star &pStar, const Model&, const std::vector<int>&, std::array<double, 2>&, std::array<double, FILTS>&);
double margEvolveWithBinary (Cluster &pCluster, Star &pStar, const Model&, const std::vector<int>&, std::array<double, 2>&, std::array<double, FILTS>&);
