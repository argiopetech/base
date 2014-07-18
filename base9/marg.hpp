#include <vector>

#include "Star.hpp"
#include "Cluster.hpp"

#include "Model.hpp"

std::vector<double> margEvolveWithBinary (const Cluster &, std::vector<StellarSystem>&, const Model&, const Isochrone&, bool);
