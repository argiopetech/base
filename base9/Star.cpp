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
    return U;
}

double Star::getMass2() const
{
    return U * massRatio;
}

void Star::setMass1(double newMass)
{
    U = newMass;
}

// *** Unused ***
// void Star::setMass2 (const Cluster &pCluster, double newMass)
// {
//     massRatio = newMass / getMass1 (pCluster);
// }

void Star::readCMD(const string &s, int filters)
{
    double tempSigma;
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

    in >> U
       >> massRatio
       >> status.at(0)
       >> clustStarPriorDens
       >> useDuringBurnIn;
}


array<double, FILTS> Star::setMags (int cmpnt, double mass, const Cluster &pCluster, const Model &evoModels, const vector<int> &filters, double &ltau)
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
        ltau = wdEvol (pCluster, evoModels, filters, mags, cmpnt);
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

double Star::wdEvol (const Cluster &pCluster, const Model &evoModels, const vector<int> &filters, array<double, FILTS> &globalMags, int cmpnt)
{
    std::pair<double, double> teffRadiusPair;

    double thisWDMass = 0.0, thisPrecLogAge = 0.0, thisLogTeff, thisWDLogRadius = 0.0;
    double mass;

    if (cmpnt == 0)
        mass = getMass1();
    else
        mass = getMass2();

    thisPrecLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(pCluster.feh, mass);

    thisWDMass = intlFinalMassReln (pCluster, evoModels, mass);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (pCluster.age, pCluster.carbonicity, thisPrecLogAge, thisWDMass);

    thisLogTeff = teffRadiusPair.first;
    thisWDLogRadius = teffRadiusPair.second;

    //*******this now gets trapped for in wdMassToTeffAndRadius so it should be unnecessary here (???)
    if (thisPrecLogAge >= pCluster.age)
    {                           // mcmc.c can cause this by adjusting masses and ages
        for (auto f : filters)
            globalMags.at(f) = -4.; // place at tip of RGB
    }
    else
    {
        //Calculate logg
        globalMags = evoModels.WDAtmosphere->teffToMags (thisLogTeff, thisWDMass, wdType.at(cmpnt));
    }

    massNow.at(cmpnt) = thisWDMass;
    wdLogTeff.at(cmpnt) = thisLogTeff;
    status.at(cmpnt) = WD;

    return thisPrecLogAge;
}

array<double, FILTS> Star::deriveCombinedMags (Matrix<double, 2, FILTS> &mag, const Cluster &pCluster, const Model &evoModels, const vector<int> &filters)
{
    auto clusterAbs = evoModels.filterSet->calcAbsCoeffs();

    assert(!filters.empty());

    double flux = 0.0;

    array<double, FILTS> combinedMags;

    combinedMags.fill(0.0);

    // can now derive combined mags
    if (mag.at(1).at(filters.front()) < 99.)
    {                           // if there is a secondary star
        for (auto f : filters)
        {
            flux = exp10((mag.at(0).at(f) / -2.5));    // add up the fluxes of the primary
            flux += exp10((mag.at(1).at(f) / -2.5));   // and the secondary
            combinedMags.at(f) = -2.5 * log10 (flux);    // (these 3 lines .at(used to?) take 5% of run time for N large)
            // if primary mag = 99.999, then this works
        }
    }                           // to make the combined mag = secondary mag
    else
    {
        for (auto f : filters)
            combinedMags.at(f) = mag.at(0).at(f);
    }

    for (decltype(filters.size()) i = 0; i < filters.size(); ++i)
    {
        int f = filters.at(i);

        combinedMags.at(f) += pCluster.mod;
        combinedMags.at(f) += (clusterAbs.at(f) - 1.0) * pCluster.abs;       // add A_.at(u-k) (standard defn of modulus already includes Av)
        photometry.at(i) = combinedMags.at(f);
    }

    return combinedMags;
}
