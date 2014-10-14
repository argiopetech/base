#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Cluster.hpp"
#include "ifmr.hpp"
#include "LinearTransform.hpp"
#include "Model.hpp"
#include "Star.hpp"

using std::ifstream;
using std::string;
using std::stringstream;
using std::vector;

vector<double> Star::msRgbEvol (const Isochrone &isochrone) const
{
    vector<double> mags;

    const auto starMass = mass;

    auto isoBegin = isochrone.eeps.cbegin();
    auto isoEnd   = isochrone.eeps.cend();

    auto m = lower_bound(isoBegin, isoEnd, starMass, EvolutionaryPoint::compareMass);

    if (m == isoEnd) {
        m -= 2;
    }
    else if (m != isoBegin) {
        m -= 1;
    }

    {
        auto magsSize = m[0].mags.size();
        mags.reserve(magsSize);

        for ( size_t f = 0; f < magsSize; ++f )
        {
            double mag = linearTransform<>(m[0].mass, m[1].mass, m[0].mags[f], m[1].mags[f], starMass).val;

            if (std::fabs(mag) < EPS)
                mags.push_back(999.99);
            else
                mags.push_back(mag);
        }
    }

    assert(mags.size() == m[0].mags.size());

    return mags;
}


vector<double> Star::getMags (const Cluster &clust, const Model &evoModels, const Isochrone &isochrone) const
{
    // Masses greater than 100Mo should never occur
    assert (mass <= 100.0);

    switch(getStatus(clust, isochrone))
    {
        case MSRG: // for main seq or giant star
            return msRgbEvol(isochrone);

        case WD:   // for white dwarf
            return wdEvol (clust, evoModels);

        default:   // For brown dwarfs, neutron stars, black hole remnants, and 0-mass secondaries
            return vector<double>(evoModels.absCoeffs.size(), 99.999);
    }
}

StarStatus Star::getStatus(const Cluster &clust, const Isochrone &isochrone) const
{
    if (mass <= 0.0001)
    {                           // for non-existent secondaries
        return DNE;
    }
    else if (mass <= isochrone.agbTipMass())
    {                           // for main seq or giant star
        return MSRG;
    }
    else if (mass <= clust.getM_wd_up())
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

vector<double> Star::wdEvol (const Cluster &clust, const Model &evoModels) const
{
    double starMass   = mass; // Cache to avoid a cache miss
    double thisWDMass = intlFinalMassReln (clust, evoModels, starMass);

    auto precLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(clust.feh, starMass, clust.yyy);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    auto teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (clust.age, clust.carbonicity, precLogAge, thisWDMass);

    double logTeff = teffRadiusPair.first;
    // double thisWDLogRadius = teffRadiusPair.second;
    
    auto mags = evoModels.WDAtmosphere->teffToMags (logTeff, thisWDMass, wdType);

    return mags;
}


// Returns actual current mass (i.e. not zams_mass)
double Star::wdMassNow(const Cluster &clust, const Model &evoModels, const Isochrone &isochrone) const
{
    if (mass <= isochrone.agbTipMass())
    {                           // for main seq or giant star
        return mass;
    }
    else if (mass <= clust.getM_wd_up())
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

    auto precLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(clust.feh, mass, clust.yyy);

    auto teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (clust.age, clust.carbonicity, precLogAge, thisWDMass);

    return teffRadiusPair.first;
}


double Star::getLtau(const Cluster &clust, const Model &evoModels) const
{
    return evoModels.mainSequenceEvol->wdPrecLogAge(clust.feh, mass, clust.yyy);
}

double StellarSystem::getMassRatio() const
{
    if (primary.mass == 0.0)
        return 0;
    else
        return secondary.mass / primary.mass;
}

void StellarSystem::setMassRatio(double r)
{
    // Secondary stars should always be between 0 and the primary's mass.
    assert ( r >= 0.0 );
    assert ( r <= 1.0 );

    secondary.mass = primary.mass * r;
}

void StellarSystem::readCMD(const string &s, int filters)
{
    int status;
    double tempSigma, massRatio;
    string starID;

    stringstream in(s);  

    in >> starID;

    for (int i = 0; i < filters; i++)
    {
        double t;
        in >> t;

        obsPhot.push_back(t);
    }

    for (int i = 0; i < filters; i++)
    {
        in >> tempSigma;

        variance.push_back(tempSigma * fabs (tempSigma));
        // The fabs() keeps the sign of the variance the same as that input by the user for sigma
        // Negative sigma (variance) is used to signal "don't count this band for this star"
    }

    in >> primary.mass
       >> massRatio
       >> status
       >> clustStarPriorDens
       >> useDuringBurnIn;

    if ((status == 1)   // MSRGB
     || (status == 3)   // WD
     || (status == 4)   // NSBH
     || (status == 5))  // BD
    {
        observedStatus = static_cast<StarStatus>(status);
    }
    else
        observedStatus = DNE;

    setMassRatio(massRatio);
}


vector<double> StellarSystem::deriveCombinedMags (const Cluster &clust, const Model &evoModels, const Isochrone &isochrone) const
{
    auto primaryMags = primary.getMags(clust, evoModels, isochrone);
    auto secondaryMags = secondary.getMags(clust, evoModels, isochrone);

    return deriveCombinedMags(clust, evoModels, isochrone, primaryMags, secondaryMags);
}

vector<double> StellarSystem::deriveCombinedMags (const Cluster &clust, const Model &evoModels, const Isochrone&, const vector<double> &primaryMags, const vector<double> &secondaryMags)
{
    assert(primaryMags.size() == secondaryMags.size());
    assert(evoModels.absCoeffs.size() == primaryMags.size());

    vector<double> combinedMags;

    auto nPrimaryMags = primaryMags.size(); // Avoid calling vector::size a lot

    // can now derive combined mags
    if (secondaryMags.front() < 99.)
    {                           // if there is a secondary star
        double flux = 0.0;

        combinedMags.reserve(nPrimaryMags);

        for (size_t f = 0; f < nPrimaryMags; ++f)
        {
            flux = exp10((primaryMags[f] / -2.5));    // add up the fluxes of the primary
            flux += exp10((secondaryMags[f] / -2.5)); // and the secondary
            combinedMags.push_back(-2.5 * __builtin_log10 (flux));    // (these 3 lines (used to?) take 5% of run time for N large)
            // if primary mag = 99.999, then this works
        }
    }  // to make the combined mag = secondary mag
    else
    {
        combinedMags = primaryMags;
    }

    for (size_t f = 0; f < nPrimaryMags; ++f)
    {
        combinedMags[f] += clust.mod;
        combinedMags[f] += (evoModels.absCoeffs[f] - 1.0) * clust.abs;       // add A_[u-k] (standard defn of modulus already includes Av)
    }

    return combinedMags;
}
