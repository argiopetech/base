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

double Star::getMass1() const
{
    return mass;
}

double Star::getMass2() const
{
    return mass2;
}

double Star::getMassRatio() const
{
    return mass2 / mass;
}

void Star::setMass1(double newMass)
{
    mass = newMass;
}

void Star::setMassRatio(double r)
{
    mass2 = mass * r;
}

// *** Unused ***
// void Star::setMass2 (const Cluster &pCluster, double newMass)
// {
//     massRatio = newMass / getMass1 (pCluster);
// }

void Star::readCMD(const string &s, int filters)
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

    in >> mass
       >> massRatio
       >> status.at(0)
       >> clustStarPriorDens
       >> useDuringBurnIn;

    setMassRatio(massRatio);
}


array<double, FILTS> Star::setMags (int cmpnt, double mass, const Cluster &pCluster, const Model &evoModels, const vector<int> &filters)
{
    array<double, FILTS> mags;

    if (mass <= 0.0001)
    {                           // for non-existent secondary stars
        for (auto f : filters)
            mags.at(f) = 99.999;

        status.at(cmpnt) = DNE;
        massNow.at(cmpnt) = 0.0;
    }
    else if (mass <= pCluster.AGBt_zmass)
    {                           // for main seq or giant star
        massNow.at(cmpnt) = mass;

        mags = evoModels.mainSequenceEvol->msRgbEvol(filters, mass);

        status.at(cmpnt) = MSRG;    // keep track of evolutionary state
    }
    else if (mass <= pCluster.M_wd_up)
    {                           // for white dwarf
        mags = wdEvol (pCluster, evoModels, cmpnt);
    }
    else if (mass <= 100.)
    {                           // for neutron star or black hole remnant
        for (auto f : filters)
            mags.at(f) = 99.999;
        status.at(cmpnt) = NSBH;
    }
    else
    {
        //     log <<  (" This condition should not happen, %.2f greater than 100 Mo\n", mass);
        for (auto f : filters)
            mags.at(f) = 99.999;

        status.at(cmpnt) = DNE;
    }

    return mags;
}

array<double, FILTS> Star::wdEvol (const Cluster &pCluster, const Model &evoModels, int cmpnt)
{
    double mass;

    if (cmpnt == 0)
        mass = getMass1();
    else
        mass = getMass2();

    double thisWDMass = intlFinalMassReln (pCluster, evoModels, mass);

    auto precLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(pCluster.feh, mass);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    auto teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (pCluster.age, pCluster.carbonicity, precLogAge, thisWDMass);

    double logTeff = teffRadiusPair.first;
    // double thisWDLogRadius = teffRadiusPair.second;
    
    auto mags = evoModels.WDAtmosphere->teffToMags (logTeff, thisWDMass, wdType.at(cmpnt));

    massNow.at(cmpnt) = thisWDMass;
    wdLogTeff.at(cmpnt) = logTeff;
    status.at(cmpnt) = WD;

    return mags;
}

double Star::getLtau(const Cluster &pCluster, const Model &evoModels, int cmpnt) const
{
    double mass;

    if (cmpnt == 0)
        mass = getMass1();
    else
        mass = getMass2();

    return evoModels.mainSequenceEvol->wdPrecLogAge(pCluster.feh, mass);
}
