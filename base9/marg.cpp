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
#include "Utility.hpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::mutex;
using std::vector;
using base::utility::ThreadPool;

// ASSUMPTIONS
//   1. margEvolveWithBinary is not adversely affected if the StellarSystem is modified
//   2. margEvolveWithBinary passes the StellarSystem with primary.mass already set appropriately
//   3. margEvolveWithBinary expects secondary.mass to be overwritten, possibly immediately
static void calcPost (const double dMass, const Cluster &clust, vector<StellarSystem> &systems, const Model &evoModels, const Isochrone &isochrone, vector<double> &post, const vector<vector<double>> &secondaryMags)
{
    // This struct is only used inside this function, so we don't want it escaping
    struct diffStruct
    {
        diffStruct(unsigned short filter, double magLower, double magUpper)
            : filter(filter), magLower(magLower), magUpper(magUpper)
        {;}
        ~diffStruct() {;}

        const unsigned short filter;
        const double magLower;
        const double magUpper;
    };

    struct starDiff
    {
        starDiff(unsigned index, vector<diffStruct> diffs)
            : index(index), diffs(diffs)
        {;}
        ~starDiff() {;}

        const unsigned index;
        const vector<diffStruct> diffs;
    };

    // Useful constants
    const auto nSystems   = systems.size();
    const auto obsSize    = systems.front().obsPhot.size();

    const double logdMass = log (dMass); // This is a pre-optimzation so we only have to call log once (instead of nEntries times)

    // Rename system.primary.mass to primaryMass for this function
    const double primaryMass = systems.front().primary.mass;
    const double dMassRatio  = log (isochrone.eeps[0].mass / primaryMass);
    const double logPrior    = clust.logPriorMass (primaryMass);

    // Get the mags based on the primary's mass set by margEvolve
    // We call deriveCombinedMags here rather than primary.getMags so that we get abs/distMod correction
    // This keeps models which haven't been converted to absolute mags from breaking
    systems.front().secondary.mass = 0.0;
    const auto primaryMags  = systems.front().primary.getMags (clust, evoModels, isochrone);

    /**** now see if any binary companions are OK ****/
    const double nSD = 3.0;           // num of st dev from obs that will contribute to likelihood

    vector<starDiff> diffs;

    // Pre-calculate diffLow and diffUp so they can be cached for the loop below
    // In its own block to contain obsFilt and obsSize
    {
        vector<diffStruct> singleDiffs; // Keeps track of diffStructs for every filter
        singleDiffs.reserve(obsSize);

        const auto combinedMags = systems.front().deriveCombinedMags (clust, evoModels, isochrone, primaryMags, secondaryMags.back());

        // expTerms is pulled out of the loop, as it's the same across every star.
        vector<double> expTerms;

        for (size_t f = 0; f < obsSize; ++f)
        {
            expTerms.push_back(exp10 ((combinedMags[f] / -2.5)));
        }

        for (size_t s = 0; s < nSystems; ++s)
        {
            singleDiffs.clear();

            auto &system = systems[s];

            for (size_t f = 0; f < obsSize; ++f)
            {
                if (system.variance[f] >= 0)
                {
                    // Optimizations(?)
                    double obsPhot = system.obsPhot[f];               // Reduces array dereference penalty
                    double nSDTerm = nSD * sqrt (system.variance[f]); // Array dereference + multiplication and sqrt call

                    const double diffLow = exp10 (((obsPhot - nSDTerm) / -2.5)) - expTerms[f];
                    const double diffUp  = exp10 (((obsPhot + nSDTerm) / -2.5)) - expTerms[f];

                    if (diffLow <= 0.0 || diffLow == diffUp) // log(-x) == NaN
                    {
                        // Instead of returning post (a half-finished calculation), return 0 to signify
                        // that no mass combinations would lead to a positive posterior density
                        singleDiffs.clear();
                        break;
                    }
                    else
                    {
                        const double magLower = -2.5 * __builtin_log10 (diffLow);
                        const double magUpper = -2.5 * __builtin_log10 (diffUp);

                        assert(!std::isnan(magLower)); // Even though we checked it above, we still check it here
                        // magUpper is allowed to be NaN, as it is checked below.

                        singleDiffs.emplace_back(f, magLower, magUpper); // Make a diffStruct and stick it in diffs
                    }
                }
            }

            // Secondary mass to 0.0 to ensure we only look at the primary for the combinedMags.
            system.primary.mass   = primaryMass;
            system.secondary.mass = 0.0;

            if ( ! singleDiffs.empty() )
            {
                double tmpLogPost = system.logPost (clust, evoModels, isochrone, logPrior, primaryMags, secondaryMags.back())
                                  + logdMass + dMassRatio;

                post[s] += __builtin_exp (tmpLogPost);

                diffs.emplace_back(s, singleDiffs);
            }
        }
    }

    // Evaluate the posterior at all reasonable points in the isochrone

    auto isoSize = isochrone.eeps.size();
    vector<double> combinedMags;

    for (size_t i = 0; i < isoSize - 1; ++i)
    {
        // More optimization. Keeps us from dereferencing this array ~4 times.
        double isoMass = isochrone.eeps[i].mass;

        if (isoMass > primaryMass)
            break;

        if (!combinedMags.empty())
            combinedMags.clear();

        for (const auto &s : diffs )
        {
            StellarSystem &system = systems[s.index];

            // All filters should be acceptable if we are to evaluate at a grid point
            // okMass starts out true, and gets set false if any filter is not acceptable
            bool okMass = true;

            // For all the diffSets (i.e., every filter)
            for (const auto &diff : s.diffs)
            {
                // The isochrone mag should be between magLower and magUpper
                // and the isochrone mass should be less than primaryMass
                //
                // Exception: If magUpper is NaN (i.e., it was < 0 as diffUp above), consider it irrelevant
                double mag = isochrone.eeps[i].mags[diff.filter];

                if ( ! (mag >= diff.magLower
                        && (mag <= diff.magUpper || std::isnan(diff.magUpper))))
                {
                    okMass = false; // If it isn't, that isn't OK
                    break;          // And break out of the loop to avoid further calculation (pre-optimization)
                }
            }

            // Now, only if we left the above loop without triggering the condition, calculate the posterior density for the system
            if (okMass)
            {
                system.setMassRatio (isoMass / primaryMass); // Set the massRatio (and therefore the secondary mass)

                if (combinedMags.empty())
                    combinedMags = system.deriveCombinedMags (clust, evoModels, isochrone, primaryMags, secondaryMags[i]);

                // Calculate the posterior density for this system
                double tmpLogPost = system.logPost (clust, evoModels, isochrone, logPrior, combinedMags)
                    + logdMass
                    + __builtin_log ((isochrone.eeps[i + 1].mass - isoMass) / primaryMass);

                post[s.index] += __builtin_exp (tmpLogPost); // And add it to the running total
            }
        }
    }
}


/* evaluate on a grid of primary mass and mass ratio to approximate the integral */
vector<double> margEvolveWithBinary (const Cluster &clust, vector<StellarSystem> &systems, const Model &evoModels, const Isochrone &isochrone, ThreadPool &pool)
{
    assert(isochrone.eeps.size() >= 2);

    vector<vector<double>> secondaryMags;

    for (const auto &e : isochrone.eeps)
    {
        Star s;
        s.mass = e.mass;

        secondaryMags.push_back(s.getMags(clust, evoModels, isochrone));
    }

    // secondaryMags.back() is the 0-mass magnitudes
    {
        Star s;
        s.mass = 0.0;

        secondaryMags.push_back(s.getMags(clust, evoModels, isochrone));
    }

    const int isoIncrem = 80;    /* ok for YY models? */

    vector<double> post(systems.size(), 0.0);

    {
        auto isoSize = isochrone.eeps.size();

        for (size_t m = 0; m < isoSize - 1; ++m)
        {
            double dIsoMass = isochrone.eeps[m + 1].mass - isochrone.eeps[m].mass;

            // In the event that we have an invalid range, skip that range
            // This generally occurs only at very high EEPs, where the masses are close together
            if (dIsoMass < 0.0)
                continue;

            double dMass    = dIsoMass / isoIncrem;

            for (auto k = 0; k < isoIncrem; ++k)
            {
                systems.front().primary.mass = isochrone.eeps[m].mass + (k * dMass);

                calcPost (dMass, clust, systems, evoModels, isochrone, post, secondaryMags);
            }
        }
    }

    return post;
}
