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

double calcPost (double, array<double, 2>&, const Cluster&, Star, const Model&, const vector<int>&, array<double, 2>&, const array<double, FILTS>&, const array<double, FILTS>&);

/* evaluate on a grid of primary mass and mass ratio to approximate the integral */
double margEvolveWithBinary (const Cluster &pCluster, const Star &pStar, const Model &evoModels, const vector<int> &filters, array<double, 2> &ltau, const array<double, FILTS> &filterPriorMin, const array<double, FILTS> &filterPriorMax)
{
    array<double, 2> mass;

    mass.at(0) = 0.0;
    mass.at(1) = 0.0;

    const struct globalIso &isochrone = evoModels.mainSequenceEvol->getIsochrone();
    double post = 0.0;

    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (pCluster.AGBt_zmass < EPS)
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
            if (dIsoMass > 0.0)
            {
                dMass = dIsoMass / isoIncrem;
                mass.at(0) = isochrone.mass.at(m) + k * dMass;

                post += calcPost (dMass, mass, pCluster, pStar, evoModels, filters, ltau, filterPriorMin, filterPriorMax);
            }
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

double calcPost (double dMass, array<double, 2> &mass, const Cluster &pCluster, Star pStar, const Model &evoModels, const vector<int> &filters, array<double, 2> &ltau, const array<double, FILTS> &filterPriorMin, const array<double, FILTS> &filterPriorMax)
{
    const struct globalIso &isochrone = evoModels.mainSequenceEvol->getIsochrone();
    double post = 0.0;

    Matrix<double, 2, FILTS> mag;

    pStar.setMass1 (mass.at(0));

    mag.at(0) = pStar.setMags (0, mass.at(0), pCluster, evoModels, filters, ltau.at(0));

    double tmpLogPost, tmpPost;

    /* first try 0.0 massRatio */
    pStar.massRatio = 0.0;

    pStar.massNow.at(1) = 0.0;
    ltau.at(1) = 0.0;          // may not be a WD, so no precursor age,
    pStar.wdLogTeff.at(1) = 0.0;      // no WD Teff,
    mag.at(1) = pStar.setMags (1, mass.at(1), pCluster, evoModels, filters, ltau.at(1));

    pStar.deriveCombinedMags (mag, pCluster, evoModels, filters);
    tmpLogPost = logPost1Star (pStar, pCluster, evoModels, filterPriorMin, filterPriorMax);
    tmpLogPost += log (dMass);
    tmpLogPost += log (isochrone.mass.at(0) / mass.at(0));    /* dMassRatio */
    tmpPost = exp (tmpLogPost);
    post += tmpPost;


    /**** now see if any binary companions are OK ****/
    double nSD = 4.0;           /* num of st dev from obs that will contribute to likelihood */
    int obsFilt = 0;

    bool isOverlap = true;          /* do the allowable masses in each filter overlap? */
    array<int, MAX_ENTRIES> okMass = { 1 };
    for (auto f : filters)
    {
        if (isOverlap)
        {
            double diffLow = exp10(((pStar.obsPhot.at(obsFilt) - nSD * sqrt (pStar.variance.at(obsFilt))) / -2.5)) - exp10((mag.at(0).at(f) / -2.5));
            double diffUp = exp10(((pStar.obsPhot.at(obsFilt) + nSD * sqrt (pStar.variance.at(obsFilt))) / -2.5)) - exp10((mag.at(0).at(f) / -2.5));
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
                    if (isochrone.mag.at(i).at(f) >= magLower && isochrone.mag.at(i).at(f) <= magUpper && isochrone.mass.at(i) <= mass.at(0))
                    {
                        okMass.at(i) *= 1; /* this mass is still ok */
                    }
                    else
                    {
                        okMass.at(i) = 0;
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
            pStar.massRatio = mass.at(0) / isochrone.mass.at(i);
            pStar.massNow.at(1) = 0.0;
            ltau.at(1) = 0.0;  // may not be a WD, so no precursor age,
            pStar.wdLogTeff.at(1) = 0.0;      // no WD Teff,
            mag.at(1) = pStar.setMags (1, mass.at(1), pCluster, evoModels, filters, ltau.at(1));

            pStar.deriveCombinedMags (mag, pCluster, evoModels, filters);
            /* now have magnitudes, want posterior probability */
            tmpLogPost = logPost1Star (pStar, pCluster, evoModels, filterPriorMin, filterPriorMax);
            tmpLogPost += log (dMass);
            tmpLogPost += log ((isochrone.mass.at(i + 1) - isochrone.mass.at(i)) / mass.at(0));
            tmpPost = exp (tmpLogPost);

            post += tmpPost;
        }
    }

    return post;
}
