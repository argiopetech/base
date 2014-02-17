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

double StellarSystem::getMassRatio() const
{
    return secondary.mass / primary.mass;
}

void StellarSystem::setMassRatio(double r)
{
    secondary.mass = primary.mass * r;
}

// *** Unused ***
// void Star::setMass2 (const Cluster &pCluster, double newMass)
// {
//     massRatio = newMass / getMass1 (pCluster);
// }

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


array<double, FILTS> Star::getMags (double mass, const Cluster &pCluster, const Model &evoModels, const vector<int> &filters)
{
    array<double, FILTS> mags;

    if (mass <= 0.0001)
    {                           // for non-existent secondary stars
        for (auto f : filters)
            mags.at(f) = 99.999;

        status = DNE;
        massNow = 0.0;
    }
    else if (mass <= pCluster.AGBt_zmass)
    {                           // for main seq or giant star
        massNow = mass;

        mags = evoModels.mainSequenceEvol->msRgbEvol(filters, mass);

        status = MSRG;    // keep track of evolutionary state
    }
    else if (mass <= pCluster.M_wd_up)
    {                           // for white dwarf
        mags = wdEvol (pCluster, evoModels);
    }
    else if (mass <= 100.)
    {                           // for neutron star or black hole remnant
        for (auto f : filters)
            mags.at(f) = 99.999;
        status = NSBH;
    }
    else
    {
        //     log <<  (" This condition should not happen, %.2f greater than 100 Mo\n", mass);
        for (auto f : filters)
            mags.at(f) = 99.999;

        status = DNE;
    }

    return mags;
}

array<double, FILTS> Star::wdEvol (const Cluster &pCluster, const Model &evoModels)
{
    double thisWDMass = intlFinalMassReln (pCluster, evoModels, mass);

    auto precLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(pCluster.feh, mass);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    auto teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (pCluster.age, pCluster.carbonicity, precLogAge, thisWDMass);

    double logTeff = teffRadiusPair.first;
    // double thisWDLogRadius = teffRadiusPair.second;
    
    auto mags = evoModels.WDAtmosphere->teffToMags (logTeff, thisWDMass, wdType);

    massNow = thisWDMass;
    wdLogTeff = logTeff;
    status = WD;

    return mags;
}

double Star::getLtau(const Cluster &pCluster, const Model &evoModels) const
{
    return evoModels.mainSequenceEvol->wdPrecLogAge(pCluster.feh, mass);
}
