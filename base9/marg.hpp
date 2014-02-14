#include <array>

#include "Star.hpp"
#include "Cluster.hpp"

#include "Model.hpp"

void deriveCombinedMags (Matrix<double, 3, FILTS> &mag, double clusterAv, double &flux, const Cluster &pCluster, Star &pStar, const Model &evoModels, const std::vector<int>&);
double margEvolveWithBinary (const Cluster &pCluster, const Star &pStar, const Model&, const std::vector<int>&, std::array<double, 2>&, std::array<double, FILTS>&, const std::array<double, FILTS>&, const std::array<double, FILTS>&);
