#include <array>
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

using std::array;
using std::vector;

const int MAX_ENTRIES = 370;

double calcPost (double, array<double, 2>&, const Cluster&, StellarSystem, const Model&, const vector<int>&);

/* evaluate on a grid of primary mass and mass ratio to approximate the integral */
double margEvolveWithBinary (const Cluster &clust, const StellarSystem &system, const Model &evoModels, const vector<int> &filters)
{
    array<double, 2> mass;

    mass.at(0) = 0.0;
    mass.at(1) = 0.0;

    const struct globalIso &isochrone = evoModels.mainSequenceEvol->getIsochrone();
    double post = 0.0;

    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (clust.AGBt_zmass < EPS)
    {
        throw WDBoundsError("Bounds error in marg.cpp");
    }

    double dMass;

    double dIsoMass = 0.0;

    int isoIncrem = 80;    /* ok for YY models? */

    assert(isochrone.nEntries >= 2);

    for (decltype(isochrone.nEntries) m = 0; m < isochrone.nEntries - 2; m++)
    {
        for (auto k = 0; k < isoIncrem; k += 1)
        {

            dIsoMass = isochrone.mass.at(m + 1) - isochrone.mass.at(m);

            /* why would dIsoMass ever be negative??? BUG in interpolation code??? */
            assert (dIsoMass > 0.0);

            dMass = dIsoMass / isoIncrem;
            mass.at(0) = isochrone.mass.at(m) + k * dMass;

            post += calcPost (dMass, mass, clust, system, evoModels, filters);
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


double calcPost (double dMass, array<double, 2> &mass, const Cluster &clust, StellarSystem system, const Model &evoModels, const vector<int> &filters)
{
    const struct globalIso &isochrone = evoModels.mainSequenceEvol->getIsochrone();
    double post = 0.0;

    array<double, FILTS> primaryMags;
    array<double, FILTS> combinedMags;

    system.primary.mass = mass.at(0);

    primaryMags = system.primary.getMags (clust, evoModels, filters);

    double tmpLogPost, tmpPost;

    /* first try 0.0 massRatio */
    system.secondary.mass = (mass.at(1));

    combinedMags = system.deriveCombinedMags (clust, evoModels, filters);
    tmpLogPost = logPost1Star (system, clust, evoModels, combinedMags);
    tmpLogPost += log (dMass);
    tmpLogPost += log (isochrone.mass.at(0) / mass.at(0));    /* dMassRatio */
    tmpPost = exp (tmpLogPost);
    post += tmpPost;

    /**** now see if any binary companions are OK ****/
    double nSD = 4.0;           /* num of st dev from obs that will contribute to likelihood */
    int obsFilt = 0;

    bool isOverlap = true;          /* do the allowable masses in each filter overlap? */
    array<bool, MAX_ENTRIES> okMass = { true };

    for (auto f : filters)
    {
        if (isOverlap)
        {
            double diffLow = exp10(((system.obsPhot.at(obsFilt) - nSD * sqrt (system.variance.at(obsFilt))) / -2.5)) - exp10((primaryMags.at(f) / -2.5));
            double diffUp = exp10(((system.obsPhot.at(obsFilt) + nSD * sqrt (system.variance.at(obsFilt))) / -2.5)) - exp10((primaryMags.at(f) / -2.5));

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
                    if ((isochrone.mag.at(i).at(f) < magLower) || (isochrone.mag.at(i).at(f) > magUpper) || (isochrone.mass.at(i) > mass.at(0)))
                    {
                        okMass.at(i) = false;
                    }
                }
            }
            obsFilt++;
        }
    }

    for (decltype(isochrone.nEntries) i = 0; i < isochrone.nEntries - 2; i++)
    {
        if (okMass.at(i))
        {
            system.setMassRatio (mass.at(0) / isochrone.mass.at(i));
            system.secondary.mass = mass.at(1);

            combinedMags = system.deriveCombinedMags (clust, evoModels, filters);

            /* now have magnitudes, want posterior probability */
            tmpLogPost = logPost1Star (system, clust, evoModels, combinedMags);
            tmpLogPost += log (dMass);
            tmpLogPost += log ((isochrone.mass.at(i + 1) - isochrone.mass.at(i)) / mass.at(0));
            tmpPost = exp (tmpLogPost);

            post += tmpPost;
        }
    }

    return post;
}
