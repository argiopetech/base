#ifndef STAR_HPP
#define STAR_HPP

#include <array>
#include <string>
#include <vector>

#include "constants.hpp"
#include "Cluster.hpp"
#include "Matrix.hpp"
#include "Model.hpp"

// Define a structure star that houses all star properties
class Star
{
  public:
    Star()
    {
        beta.fill(0.0);
    }

    ~Star() {;}

    double getLtau(const Cluster&, const Model&) const;
    double wdLogTeff(const Cluster&, const Model&) const;
    double wdMassNow(double, const Cluster&, const Model&) const;

    std::array<double, FILTS> getMags (double, const Cluster&, const Model&, const std::vector<int>&);
    std::array<double, FILTS> wdEvol (const Cluster&, const Model&) const;

    std::array<double, NPARAMS> beta;

    int status = 0;

    WdAtmosphere wdType = WdAtmosphere::DA;

    double mass = 0.0;
};

class StellarSystem
{
  public:
    StellarSystem()
    {
        obsPhot.fill(0.0);
        variance.fill(0.0);
    }

    StellarSystem(const std::string &s, int filters)
        : StellarSystem()
    {
        readCMD(s, filters);
    }

    void readCMD(const std::string&, int);

    double getMassRatio() const;
    void setMassRatio(double);

    Star primary;
    Star secondary;

    std::array<double, FILTS> obsPhot;
    std::array<double, FILTS> variance;

    bool useDuringBurnIn = false;       // switch whether to use star to burn in cluster parameters

    double clustStarPriorDens = 0.0;    // prior probability that the star is a cluster star
    double clustStarProposalDens = 0.0; // proposal density for steps to the cluster star model
};

#endif
