#ifndef STAR_HPP
#define STAR_HPP

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
    StarStatus getStatus(const Cluster&, const Isochrone&) const;

    double getLtau(const Cluster&, const Model&) const;
    double wdLogTeff(const Cluster&, const Model&) const;
    double wdMassNow(const Cluster&, const Model&, const Isochrone&) const;

    std::vector<double> getMags (const Cluster&, const Model&, const Isochrone&) const;
    std::vector<double> wdEvol (const Cluster&, const Model&) const;
    std::vector<double> msRgbEvol (const Isochrone&) const;

    // Variables
    WdAtmosphere wdType = WdAtmosphere::DA;

    double mass = 0.0;
};

class StellarSystem
{
  public:
    StellarSystem()
    {;}

    StellarSystem(const std::string &s, int filters)
        : StellarSystem()
    {
        readCMD(s, filters);
    }

    ~StellarSystem() {;}

    // Functions
    void readCMD(const std::string&, int);

    void setMassRatio(double);
    double getMassRatio() const;
    double logPost (const Cluster &clust, const Model &evoModels, const Isochrone&) const;
    double logPost (const Cluster &clust, const Model &evoModels, const Isochrone&, double, const std::vector<double>&) const;
    double logPost (const Cluster &clust, const Model &evoModels, const Isochrone&, double, const std::vector<double>&, const std::vector<double>&) const;

    std::vector<double> deriveCombinedMags (const Cluster&, const Model&, const Isochrone&) const;
    static std::vector<double> deriveCombinedMags (const Cluster&, const Model&, const Isochrone&, const std::vector<double>&, const std::vector<double>&);

    // Variables
    Star primary;
    Star secondary;

    StarStatus observedStatus;

    std::vector<double> obsPhot;
    std::vector<double> variance;

    bool useDuringBurnIn = false;       // switch whether to use star to burn in cluster parameters

    double clustStarPriorDens = 0.0;    // prior probability that the star is a cluster star
    double clustStarProposalDens = 0.0; // proposal density for steps to the cluster star model

  private:
    double lastPrimary = std::numeric_limits<double>::quiet_NaN();

    std::vector<double> lastPrimaryMags;
};


class InvalidStar : public std::range_error
{
  public:
    explicit InvalidStar (const std::string& what_arg)
        : std::range_error(what_arg) {}

    explicit InvalidStar (const char* what_arg)
        : std::range_error(what_arg) {}
};

#endif
