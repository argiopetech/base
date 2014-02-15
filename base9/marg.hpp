#include <array>

#include "Star.hpp"
#include "Cluster.hpp"

#include "Model.hpp"

double margEvolveWithBinary (const Cluster &, const Star &, const Model&, const std::vector<int>&, const std::array<double, FILTS>&, const std::array<double, FILTS>&);
std::array<double, FILTS> deriveCombinedMags (Matrix<double, 2, FILTS>&, const Cluster&, const Model&, const std::vector<int>&);
