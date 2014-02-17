#include <array>

#include "Star.hpp"
#include "Cluster.hpp"

#include "Model.hpp"

double margEvolveWithBinary (const Cluster &, const StellarSystem &, const Model&, const std::vector<int>&);
std::array<double, FILTS> deriveCombinedMags (const std::array<double, FILTS>&, const std::array<double, FILTS>&, const Cluster&, const Model&, const std::vector<int>&);
