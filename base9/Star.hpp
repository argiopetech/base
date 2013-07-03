#include <array>

#include "constants.hpp"
#include "Matrix.hpp"

// Define a structure star that houses all star properties
class Star
{
  public:
    Star()
    {
        // NPARAMS-element arrays
        useFilt.fill(false);
        obsPhot.fill(0.0);
        photometry.fill(0.0);
        variance.fill(0.0);

        // 2-element arrays
        status.fill(0);
        wdType.fill(0);
        massNow.fill(0.0);
        wdLogTeff.fill(0.0);
        betaMassRatio.fill(0.0);

        for (auto &a : beta)
            a.fill(0.0);
    }
    ~Star() {;}

    Matrix<double, NPARAMS, 2> beta;

    std::array<bool, FILTS> useFilt;

    std::array<int, 2> status;
    std::array<int, 2> wdType;

    std::array<double, FILTS> obsPhot;
    std::array<double, FILTS> photometry;
    std::array<double, FILTS> variance;
    std::array<double, 2> betaMassRatio;
    std::array<double, 2> wdLogTeff;
    std::array<double, 2> massNow;      // Actual current masses of each component (i.e. not zams_mass)

    bool isFieldStar = false;
    bool useDuringBurnIn = false;       // switch whether to use star to burn in cluster parameters

    double clustStarPriorDens = 0.0;    // prior probability that the star is a cluster star
    double clustStarProposalDens = 0.0; // proposal density for steps to the cluster star model
    double U = 0.0;
    double massRatio = 0.0;             // massRatio = secondary mass / primary mass (between 0 and 1)
    double meanU = 0.0;
    double varU = 0.0;
    double meanMassRatio = 0.0;
    double varMassRatio = 0.0;
    double UStepSize = 0.0;
    double massRatioStepSize = 0.0;
};
