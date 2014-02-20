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

#include <iostream>

using std::array;
using std::vector;

int okMasses;

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
                      + log (isochrone.mass.at(0) / primaryMass);   // dMassRatio
    double post = exp (tmpLogPost);

    /**** now see if any binary companions are OK ****/
    const double nSD = 4.0;           /* num of st dev from obs that will contribute to likelihood */

    boost::dynamic_bitset<> okMass(isochrone.nEntries, false);

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

    vector<diffStruct> diffs;

    {
        int obsFilt = 0;
    
        for (auto f : filters)
        {
            const double diffLow = exp10 (((system.obsPhot.at(obsFilt) - nSD * sqrt (system.variance.at(obsFilt))) / -2.5))
                                 - exp10 ((primaryMags.at(f) / -2.5));
            const double diffUp = exp10 (((system.obsPhot.at(obsFilt) + nSD * sqrt (system.variance.at(obsFilt))) / -2.5))
                                - exp10 ((primaryMags.at(f) / -2.5));

            if (diffLow <= 0.0 || diffUp <= 0.0 || diffLow == diffUp)
            {
                return 0.0;
            }
            else
            {
                const double magLower = -2.5 * log10 (diffLow);
                const double magUpper = -2.5 * log10 (diffUp);

                assert(!std::isnan(magLower));
                assert(!std::isnan(magUpper));
            
                diffs.emplace_back(obsFilt, magLower, magUpper);
            }

            ++obsFilt;
        }
    }

    for (decltype(isochrone.nEntries) i = 0; i < isochrone.nEntries - 1; ++i)
    {
        bool okMass = true;

        for (auto diff : diffs)
        {
            if (! (isochrone.mag.at(i).at(diff.filter) >= diff.magLower
                && isochrone.mag.at(i).at(diff.filter) <= diff.magUpper
                && isochrone.mass.at(i) <= primaryMass))
            {
                okMass = false;
                break;
            }
        }

        if (okMass)
        {
            ++okMasses;

            system.setMassRatio (isochrone.mass.at(i) / primaryMass);

            /* now have magnitudes, want posterior probability */
            tmpLogPost = system.logPost (clust, evoModels, filters)
                       + logdMass
                       + log ((isochrone.mass.at(i + 1) - isochrone.mass.at(i)) / primaryMass);

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

    const int isoIncrem = 80;    /* ok for YY models? */

    double post = 0.0;
    okMasses = 0;

    const struct globalIso &isochrone = evoModels.mainSequenceEvol->getIsochrone();
    assert(isochrone.nEntries >= 2);

    for (decltype(isochrone.nEntries) m = 0; m < isochrone.nEntries - 1; m++)
    {
        double dIsoMass = isochrone.mass.at(m + 1) - isochrone.mass.at(m);

        /* why would dIsoMass ever be negative??? BUG in interpolation code??? */
        assert (dIsoMass >= 0.0);

        double dMass = dIsoMass / isoIncrem;

        for (auto k = 0; k < isoIncrem; k += 1)
        {
            system.primary.mass = isochrone.mass.at(m) + (k * dMass);

            post += calcPost (dMass, clust, system, evoModels, filters);
        }
    }

    std::cout << okMasses << " ok Masses" << std::endl;

    if (post >= 0.0)
    {
        return post;
    }
    else
    {
        return 0.0;
    }
}
