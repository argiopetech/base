#include <vector>

#include "Star.hpp"
#include "Cluster.hpp"

#include "Model.hpp"
#include "Utility.hpp"

std::vector<double> margEvolveWithBinary (const Cluster &, std::vector<StellarSystem>&, const Model&, const Isochrone&, base::utility::ThreadPool&, bool);
