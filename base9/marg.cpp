#include <array>
#include <bitset>
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
#include "MsFilterSet.hpp"
#include "WhiteDwarf.hpp"

#include "boost/dynamic_bitset.hpp"

using std::array;
using std::vector;

const int isoIncrem = 80;    /* ok for YY models? */

// ASSUMPTIONS
//   1. margEvolveWithBinary is not adversely affected if the StellarSystem is modified
//   2. margEvolveWithBinary passes the StellarSystem with primary.mass already set appropriately
//   3. margEvolveWithBinary expects secondary.mass to be overwritten, possibly immediately
static double calcPost (const double dMass, const Cluster &clust, StellarSystem &system, const Model &evoModels, const vector<int> &filters)
{
    const struct globalIso &isochrone = evoModels.mainSequenceEvol->getIsochrone();

    // Get the mags based on the primary's mass set by margEvolve
    // Also, go ahead and rename system.primary.mass to primaryMass for this function
    array<double, FILTS> primaryMags = system.primary.getMags (clust, evoModels, filters);
    const double primaryMass = system.primary.mass;

    // Other useful constants
    const double logdMass = log (dMass);

    /* first try 0.0 massRatio */
    system.secondary.mass = 0.0;

    double tmpLogPost = system.logPost (clust, evoModels, filters)
                      + logdMass
                      + log (isochrone.mass[0] / primaryMass);   // dMassRatio
    double post = exp (tmpLogPost);

    /**** now see if any binary companions are OK ****/
    const double nSD = 4.0;           /* num of st dev from obs that will contribute to likelihood */
    int obsFilt = 0;

    bool isOverlap = true;          /* do the allowable masses in each filter overlap? */
    boost::dynamic_bitset<> okMass(isochrone.nEntries, true);

    for (auto f : filters)
    {
        if (isOverlap)
        {
            double diffLow = exp10 (((system.obsPhot[obsFilt] - nSD * sqrt (system.variance[obsFilt])) / -2.5))
                           - exp10 ((primaryMags[f] / -2.5));
            double diffUp = exp10 (((system.obsPhot[obsFilt] + nSD * sqrt (system.variance[obsFilt])) / -2.5))
                          - exp10 ((primaryMags[f] / -2.5));

            if (diffLow <= 0.0 || diffUp <= 0.0 || diffLow == diffUp)
            {
                isOverlap = false;
            }
            else
            {
                double magLower = -2.5 * log10 (diffLow);
                double magUpper = -2.5 * log10 (diffUp);

                for (decltype(isochrone.nEntries) i = 0; i < isochrone.nEntries - 1; i++)
                {
                    if (!(isochrone.mag[i][f] >= magLower && isochrone.mag[i][f] <= magUpper && isochrone.mass[i] <= primaryMass))
                    {
                        okMass[i] = false;
                    }
                }
            }
            obsFilt++;
        }
    }

    for (decltype(isochrone.nEntries) i = 0; i < isochrone.nEntries - 2; ++i)
    {
        if (okMass[i])
        {
            system.setMassRatio (isochrone.mass[i] / primaryMass);

            /* now have magnitudes, want posterior probability */
            tmpLogPost = system.logPost (clust, evoModels, filters)
                       + logdMass
                       + log ((isochrone.mass[i + 1] - isochrone.mass[i]) / primaryMass);

            post += exp (tmpLogPost);
        }
    }

    return post;
}


/* evaluate on a grid of primary mass and mass ratio to approximate the integral */
double margEvolveWithBinary (const Cluster &clust, StellarSystem system, const Model &evoModels, const vector<int> &filters)
{
    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (clust.AGBt_zmass < EPS)
    {
        throw WDBoundsError("Bounds error in marg.cpp");
    }

    double post = 0.0;

    const struct globalIso &isochrone = evoModels.mainSequenceEvol->getIsochrone();
    assert(isochrone.nEntries >= 2);

    for (decltype(isochrone.nEntries) m = 0; m < isochrone.nEntries - 2; m++)
    {
        for (auto k = 0; k < isoIncrem; k += 1)
        {

            double dIsoMass = isochrone.mass.at(m + 1) - isochrone.mass.at(m);

            /* why would dIsoMass ever be negative??? BUG in interpolation code??? */
            assert (dIsoMass >= 0.0);

            double dMass = dIsoMass / isoIncrem;
            system.primary.mass = isochrone.mass.at(m) + k * dMass;

            post += calcPost (dMass, clust, system, evoModels, filters);
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
