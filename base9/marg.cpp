#include <array>
#include <bitset>
#include <iostream>
#include <vector>

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "Cluster.hpp"
#include "Star.hpp"

#include "densities.hpp"
#include "Matrix.hpp"
#include "Model.hpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::vector;

// ASSUMPTIONS
//   1. margEvolveWithBinary is not adversely affected if the StellarSystem is modified
//   2. margEvolveWithBinary passes the StellarSystem with primary.mass already set appropriately
//   3. margEvolveWithBinary expects secondary.mass to be overwritten, possibly immediately
static double calcPost (const double dMass, const Cluster &clust, StellarSystem &system, const Model &evoModels)
{
    // This struct is only used inside this function, so we don't want it escaping
    struct diffStruct
    {
        diffStruct(int filter, double magLower, double magUpper)
            : filter(filter), magLower(magLower), magUpper(magUpper)
        {;}
        ~diffStruct() {;}

        const int filter;
        const double magLower;
        const double magUpper;
    };

    auto &isochrone = evoModels.mainSequenceEvol->getIsochrone();

    // Secondary mass to 0.0 to ensure we only look at the primary for the primaryMags.
    system.secondary.mass = 0.0;

    // Get the mags based on the primary's mass set by margEvolve
    // Also, go ahead and rename system.primary.mass to primaryMass for this function
    // We call deriveCombinedMags here rather than primary.getMags so that we get abs/distMod correction
    // This keeps models which haven't been converted to absolute mags from breaking
    vector<double> primaryMags = system.deriveCombinedMags (clust, evoModels);
    const double primaryMass = system.primary.mass;

    // Other useful constants
    const double logdMass = log (dMass); // This is a pre-optimzation so we only have to call log once (instead of nEntries times)

    double tmpLogPost = system.logPost (clust, evoModels)
                      + logdMass
                      + log (isochrone.eeps.at(0).mass / primaryMass);   // dMassRatio
    double post = exp (tmpLogPost);

    /**** now see if any binary companions are OK ****/
    const double nSD = 3.0;           // num of st dev from obs that will contribute to likelihood

    vector<diffStruct> diffs; // Keeps track of diffStructs for every filter

    // Pre-calculate diffLow and diffUp so they can be cached for the loop below
    // In its own block to contain obsFilt
    for (size_t f = 0; f < system.obsPhot.size(); ++f)
    {
        if (system.variance.at(f) >= 0)
        {
            const double diffLow = exp10 (((system.obsPhot.at(f) - nSD * sqrt (system.variance.at(f))) / -2.5))
                                 - exp10 ((primaryMags.at(f) / -2.5));
            const double diffUp = exp10 (((system.obsPhot.at(f) + nSD * sqrt (system.variance.at(f))) / -2.5))
                                - exp10 ((primaryMags.at(f) / -2.5));

            if (diffLow <= 0.0 || diffLow == diffUp) // log(-x) == NaN
            {
                return 0.0; // Instead of returning post (a half-finished calculation), return 0 to signify
                // that no mass combinations would lead to a positive posterior density
            }
            else
            {
                const double magLower = -2.5 * log10 (diffLow);
                const double magUpper = -2.5 * log10 (diffUp);

                assert(!std::isnan(magLower)); // Even though we checked it above, we still check it here
                // magUpper is allowed to be NaN, as it is checked below.

                diffs.emplace_back(f, magLower, magUpper); // Make a diffStruct and stick it in diffs
            }
        }
    }

    // This only happens if you have all variances negative for a star
    assert (!diffs.empty());

    // Evaluate the posterior at all reasonable points in the isochrone
    for (decltype(isochrone.eeps.size()) i = 0; i < isochrone.eeps.size() - 1; ++i)
    {
        // All filters should be acceptable if we are to evaluate at a grid point
        // okMass starts out true, and gets set false if any filter is not acceptable
        bool okMass = true;

        // For all the diffSets (i.e., every filter)
        for (auto diff : diffs)
        {
            // The isochrone mag should be between magLower and magUpper
            // and the isochrone mass should be less than primaryMass
            //
            // Exception: If magUpper is NaN (i.e., it was < 0 as diffUp above), consider it irrelevant
            if ( ! (isochrone.eeps.at(i).mags.at(diff.filter) >= diff.magLower
                && (isochrone.eeps.at(i).mags.at(diff.filter) <= diff.magUpper || std::isnan(diff.magUpper))
                &&  isochrone.eeps.at(i).mass <= primaryMass))
            {
                okMass = false; // If it isn't, that isn't OK
                break;          // And break out of the loop to avoid further calculation (pre-optimization)
            }
        }

        // Now, only if we left the above loop without triggering the condition, calculate the posterior density for the system
        if (okMass)
        {
            system.setMassRatio (isochrone.eeps.at(i).mass / primaryMass); // Set the massRatio (and therefore the secondary mass)

            // Calculate the posterior density for this system
            tmpLogPost = system.logPost (clust, evoModels)
                       + logdMass
                       + log ((isochrone.eeps.at(i + 1).mass - isochrone.eeps.at(i).mass) / primaryMass);

            post += exp (tmpLogPost); // And add it to the running total
        }
    }

    return post;
}


/* evaluate on a grid of primary mass and mass ratio to approximate the integral */
double margEvolveWithBinary (const Cluster &clust, StellarSystem system, const Model &evoModels)
{
    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (clust.AGBt_zmass < EPS)
    {
        throw WDBoundsError("Bounds error in marg.cpp");
    }

    const int isoIncrem = 80;    /* ok for YY models? */

    double post = 0.0;

    auto &isochrone = evoModels.mainSequenceEvol->getIsochrone();
    assert(isochrone.eeps.size() >= 2);

    for (decltype(isochrone.eeps.size()) m = 0; m < isochrone.eeps.size() - 1; m++)
    {
        double dIsoMass = isochrone.eeps.at(m + 1).mass - isochrone.eeps.at(m).mass;

        /* why would dIsoMass ever be negative??? BUG in interpolation code??? */
        assert (dIsoMass >= 0.0);

        double dMass = dIsoMass / isoIncrem;

        for (auto k = 0; k < isoIncrem; k += 1)
        {
            system.primary.mass = isochrone.eeps.at(m).mass + (k * dMass);

            post += calcPost (dMass, clust, system, evoModels);
        }
    }

    if (post >= 0.0)
    {
        return post;
    }
    else
    {
        return 0.0;
    }
}
