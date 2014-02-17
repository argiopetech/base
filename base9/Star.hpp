#ifndef STAR_HPP
#define STAR_HPP

#include <array>
#include <string>
#include <vector>

#include "constants.hpp"
#include "Cluster.hpp"
#include "Matrix.hpp"
#include "Model.hpp"

class Star
{
  public:
    Star() {;}
    ~Star() {;}

    // Functions
    int getStatus(const Cluster&) const;

    double getLtau(const Cluster&, const Model&) const;
    double wdLogTeff(const Cluster&, const Model&) const;
    double wdMassNow(const Cluster&, const Model&) const;

    std::array<double, FILTS> getMags (const Cluster&, const Model&, const std::vector<int>&) const;
    std::array<double, FILTS> wdEvol (const Cluster&, const Model&) const;


    // Variables
    WdAtmosphere wdType = WdAtmosphere::DA;

    int status = 0;

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

    ~StellarSystem() {;}


    // Functions
    void readCMD(const std::string&, int);

    double getMassRatio() const;
    void setMassRatio(double);


    // Variables
    Star primary;
    Star secondary;

    std::array<double, FILTS> obsPhot;
    std::array<double, FILTS> variance;

    bool useDuringBurnIn = false;       // switch whether to use star to burn in cluster parameters

    double clustStarPriorDens = 0.0;    // prior probability that the star is a cluster star
    double clustStarProposalDens = 0.0; // proposal density for steps to the cluster star model
};

#endif
