#include <array>
#include <vector>

#include "Cluster.hpp"
#include "Star.hpp"

class Chain
{
  public:
    std::array<bool, NPARAMS> acceptClust;
    std::array<bool, NPARAMS> rejectClust;
    std::vector<StellarSystem> stars;

    Cluster clust;

    int *acceptMass;
    int *rejectMass;
    int *acceptMassRatio;
    int *rejectMassRatio;
    int *isFieldStarNow;
    int *isClusterStarNow;

    double temperature;
};
