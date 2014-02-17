#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "Cluster.hpp"
#include "ifmr.hpp"
#include "Model.hpp"
#include "Star.hpp"

using std::array;
using std::ifstream;
using std::string;
using std::stringstream;
using std::vector;

array<double, FILTS> Star::getMags (const Cluster &pCluster, const Model &evoModels, const vector<int> &filters) const
{
    array<double, FILTS> mags;

    if (mass <= 0.0001)
    {                           // for non-existent secondary stars
        for (auto f : filters)
            mags.at(f) = 99.999;
    }
    else if (mass <= pCluster.AGBt_zmass)
    {                           // for main seq or giant star
        mags = evoModels.mainSequenceEvol->msRgbEvol(filters, mass);
    }
    else if (mass <= pCluster.M_wd_up)
    {                           // for white dwarf
        mags = wdEvol (pCluster, evoModels);
    }
    else if (mass <= 100.)
    {                           // for neutron star or black hole remnant
        for (auto f : filters)
            mags.at(f) = 99.999;
    }
    else
    {
        //     log <<  (" This condition should not happen, %.2f greater than 100 Mo\n", mass);
        for (auto f : filters)
            mags.at(f) = 99.999;
    }

    return mags;
}

int Star::getStatus(const Cluster &clust) const
{
    if (mass <= clust.AGBt_zmass)
    {                           // for main seq or giant star
        return MSRG;
    }
    else if (mass <= clust.M_wd_up)
    {                           // for white dwarf
        return WD;
    }
    else if (mass <= 100.)
    {                           // for neutron star or black hole remnant
        return NSBH;
    }
    else
    {                           // for everything else
        return DNE;
    }
}

array<double, FILTS> Star::wdEvol (const Cluster &clust, const Model &evoModels) const
{
    double thisWDMass = intlFinalMassReln (clust, evoModels, mass);

    auto precLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(clust.feh, mass);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    auto teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (clust.age, clust.carbonicity, precLogAge, thisWDMass);

    double logTeff = teffRadiusPair.first;
    // double thisWDLogRadius = teffRadiusPair.second;
    
    auto mags = evoModels.WDAtmosphere->teffToMags (logTeff, thisWDMass, wdType);

    return mags;
}


// Returns actual current mass (i.e. not zams_mass)
double Star::wdMassNow(const Cluster &clust, const Model &evoModels) const
{
    if (mass <= clust.AGBt_zmass)
    {                           // for main seq or giant star
        return mass;
    }
    else if (mass <= clust.M_wd_up)
    {                           // for white dwarf
        return intlFinalMassReln (clust, evoModels, mass);
    }
    else
    {
        return 0.0;
    }
}


double Star::wdLogTeff(const Cluster &clust, const Model &evoModels) const
{
    double thisWDMass = intlFinalMassReln (clust, evoModels, mass);

    auto precLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(clust.feh, mass);

    auto teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (clust.age, clust.carbonicity, precLogAge, thisWDMass);

    return teffRadiusPair.first;
}


double Star::getLtau(const Cluster &clust, const Model &evoModels) const
{
    return evoModels.mainSequenceEvol->wdPrecLogAge(clust.feh, mass);
}

double StellarSystem::getMassRatio() const
{
    return secondary.mass / primary.mass;
}

void StellarSystem::setMassRatio(double r)
{
    secondary.mass = primary.mass * r;
}

void StellarSystem::readCMD(const string &s, int filters)
{
    double tempSigma, massRatio;
    string starID;

    stringstream in(s);  

    in >> starID;

    for (int i = 0; i < filters; i++)
    {
        in >> obsPhot.at(i);
    }

    for (int i = 0; i < filters; i++)
    {
        in >> tempSigma;

        variance.at(i) = tempSigma * fabs (tempSigma);
        // The fabs() keeps the sign of the variance the same as that input by the user for sigma
        // Negative sigma (variance) is used to signal "don't count this band for this star"
    }

    in >> primary.mass
       >> massRatio
       >> primary.status
       >> clustStarPriorDens
       >> useDuringBurnIn;

    setMassRatio(massRatio);
}
