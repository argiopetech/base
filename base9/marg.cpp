#include <array>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "evolve.hpp"
#include "gBaraffeMag.hpp"
#include "densities.hpp"
#include "Model.hpp"
#include "wdEvol.hpp"
#include "MsFilterSet.hpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::vector;

const int MAX_ENTRIES = 370;

void calcPost (double*, double, double[][FILTS], double, double*, double*, const Cluster&, Star&, const Model&, const vector<int>&, array<double, 2>&, array<double, FILTS>&, const array<double, FILTS>&, const array<double, FILTS>&);

/* evaluate on a grid of primary mass and mass ratio to approximate the integral */
double margEvolveWithBinary (const Cluster &pCluster, const Star &pStar, const Model &evoModels, const vector<int> &filters, array<double, 2> &ltau, array<double, FILTS> &globalMags, const array<double, FILTS> &filterPriorMin, const array<double, FILTS> &filterPriorMax)
{
    double mag[3][FILTS], mass[2], flux, clusterAv;
    double post = 0.0;

    mass[0] = 0.0;
    mass[1] = 0.0;

    const struct globalIso &isochrone = evoModels.mainSequenceEvol->getIsochrone();

    //Don't recalculate AGB mass (and isochrone) if these parameters are the same as they
    //were last time through
    // if (fabs (isochrone.FeH - pCluster.getFeH()) > EPS || fabs (isochrone.logAge - pCluster.getAge()) > EPS || fabs (isochrone.Y - pCluster.getY()) > EPS)
    // {
    //     pCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(filters, pCluster.getFeH(), pCluster.getY(), pCluster.getAge());        // determine AGBt ZAMS mass, to find evol state
    // }

    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (pCluster.AGBt_zmass < EPS)
    {
        throw WDBoundsError("Bounds error in marg.cpp");
    }

    clusterAv = pCluster.getAbs();

    double dMass;

    double dIsoMass = 0.0;

    double isoIncrem = 80.0;    /* ok for YY models? */

    Star myStar(pStar);

    for (int m = 0; m < isochrone.nEntries - 2; m++)
    {
        for (int k = 0; k < isoIncrem; k += 1)
        {

            dIsoMass = isochrone.mass[m + 1] - isochrone.mass[m];

            /* why would dIsoMass ever be negative??? BUG in interpolation code??? */
            if (dIsoMass > 0.0)
            {
                dMass = dIsoMass / isoIncrem;
                mass[0] = isochrone.mass[m] + k * dMass;

                calcPost (&post, dMass, mag, clusterAv, &flux, mass, pCluster, myStar, evoModels, filters, ltau, globalMags, filterPriorMin, filterPriorMax);
            }
        }
    }
    if (post > 0.0)
    {
        return post;
    }
    else
    {
        return 0.0;
    }
}

void setMags (double mag[][FILTS], int cmpnt, double *mass, const Cluster &pCluster, Star &pStar, const Model &evoModels, const vector<int> &filters,  array<double, 2> &ltau, array<double, FILTS> &globalMags)
{
    if (mass[cmpnt] <= 0.0001)
    {                           // for non-existent secondary stars
        for (auto f : filters)
            mag[cmpnt][f] = 99.999;
        pStar.status[cmpnt] = DNE;
        pStar.massNow[cmpnt] = 0.0;
    }
    else if (mass[cmpnt] <= pCluster.AGBt_zmass)
    {                           // for main seq or giant star
        pStar.massNow[cmpnt] = evoModels.mainSequenceEvol->msRgbEvol(filters, globalMags, mass[cmpnt]);
        for (auto f : filters)
            mag[cmpnt][f] = globalMags[f];
        pStar.status[cmpnt] = MSRG;    // keep track of evolutionary state
    }
    else if (mass[cmpnt] <= pCluster.M_wd_up)
    {                           // for white dwarf
        ltau[cmpnt] = wdEvol (pCluster, evoModels, filters, globalMags, pStar, cmpnt);
        for (auto f : filters)
            mag[cmpnt][f] = globalMags[f];
    }
    else if (mass[cmpnt] <= 100.)
    {                           // for neutron star or black hole remnant
        for (auto f : filters)
            mag[cmpnt][f] = 99.999;
        pStar.status[cmpnt] = NSBH;
    }
    else
    {
        //     log <<  (" This condition should not happen, %.2f greater than 100 Mo\n", mass[cmpnt]);
        for (auto f : filters)
                mag[cmpnt][f] = 99.999;
        pStar.status[cmpnt] = DNE;
    }
}

void deriveCombinedMags (double mag[][FILTS], double clusterAv, double *flux, const Cluster &pCluster, Star &pStar, const Model &evoModels, const vector<int> &filters)
{
    auto clusterAbs = evoModels.filterSet->calcAbsCoeffs();

    // can now derive combined mags
    if (mag[1][filters.front()] < 99.)
    {                           // if there is a secondary star
        for (auto f : filters)
        {
                (*flux) = pow (10.0, (mag[0][f] / -2.5));    // add up the fluxes of the primary
                (*flux) += pow (10.0, (mag[1][f] / -2.5));   // and the secondary
                mag[2][f] = -2.5 * log10 (*flux);    // (these 3 lines take 5% of run time for N large)
                // if primary mag = 99.999, then this works
        }
    }                           // to make the combined mag = secondary mag
    else
    {
        for (auto f : filters)
            mag[2][f] = mag[0][f];
    }

    for (auto f : filters)
    {
            mag[2][f] += pCluster.getMod();
            mag[2][f] += (clusterAbs[f] - 1.0) * clusterAv;       // add A_[u-k] (standard defn of modulus already includes Av)
            pStar.photometry[f] = mag[2][f];
    }
}


void calcPost (double *post, double dMass, double mag[][FILTS], double clusterAv, double *flux, double *mass, const Cluster &pCluster, Star &pStar, const Model &evoModels, const vector<int> &filters, array<double, 2> &ltau, array<double, FILTS> &globalMags, const array<double, FILTS> &filterPriorMin, const array<double, FILTS> &filterPriorMax)
{
    const struct globalIso &isochrone = evoModels.mainSequenceEvol->getIsochrone();

    pStar.setMass1 (pCluster, mass[0]);

    int cmpnt = 0;

    setMags (mag, cmpnt, mass, pCluster, pStar, evoModels, filters, ltau, globalMags);

    double tmpLogPost, tmpPost;

    /* first try 0.0 massRatio */
    cmpnt = 1;
    pStar.massRatio = 0.0;

    for (auto f : filters)
            globalMags[f] = 99.999;

    pStar.massNow[cmpnt] = 0.0;
    ltau[cmpnt] = 0.0;          // may not be a WD, so no precursor age,
    pStar.wdLogTeff[cmpnt] = 0.0;      // no WD Teff,
    setMags (mag, cmpnt, mass, pCluster, pStar, evoModels, filters, ltau, globalMags);

    deriveCombinedMags (mag, clusterAv, flux, pCluster, pStar, evoModels, filters);
    tmpLogPost = logPost1Star (pStar, pCluster, evoModels, filterPriorMin, filterPriorMax);
    tmpLogPost += log (dMass);
    tmpLogPost += log (isochrone.mass[0] / mass[0]);    /* dMassRatio */
    tmpPost = exp (tmpLogPost);
    (*post) += tmpPost;


    /**** now see if any binary companions are OK ****/
    double magLower;
    double magUpper;
    double nSD = 4.0;           /* num of st dev from obs that will contribute to likelihood */
    int obsFilt = 0;
    int i;

    bool isOverlap = true;          /* do the allowable masses in each filter overlap? */
    double diffLow = 0.0;
    double diffUp = 0.0;
    int okMass[MAX_ENTRIES] = { 1 };
    for (auto f : filters)
    {
        if (isOverlap)
        {
            diffLow = pow (10.0, ((pStar.obsPhot[obsFilt] - nSD * sqrt (pStar.variance[obsFilt])) / -2.5)) - pow (10.0, (mag[0][f] / -2.5));
            diffUp = pow (10.0, ((pStar.obsPhot[obsFilt] + nSD * sqrt (pStar.variance[obsFilt])) / -2.5)) - pow (10.0, (mag[0][f] / -2.5));
            if (diffLow <= 0.0 || diffUp <= 0.0 || diffLow == diffUp)
            {
                isOverlap = false;

                /**** necessary here??? *****/
                for ( auto f : filters )
                    mag[2][f] = mag[0][f];
            }
            else
            {
                magLower = -2.5 * log10 (diffLow);
                magUpper = -2.5 * log10 (diffUp);

                for (i = 0; i < isochrone.nEntries - 1; i++)
                {
                    if (isochrone.mag[i][f] >= magLower && isochrone.mag[i][f] <= magUpper && isochrone.mass[i] <= mass[0])
                    {
                        okMass[i] *= 1; /* this mass is still ok */
                    }
                    else
                    {
                        okMass[i] = 0;
                    }
                }
            }
            obsFilt++;
        }
    }

    for (i = 0; i < isochrone.nEntries - 2; i++)
    {
        if (okMass[i])
        {
            cmpnt = 1;
            pStar.massRatio = mass[0] / isochrone.mass[i];
            for (auto f : filters)
                globalMags[f] = 99.999;
            pStar.massNow[cmpnt] = 0.0;
            ltau[cmpnt] = 0.0;  // may not be a WD, so no precursor age,
            pStar.wdLogTeff[cmpnt] = 0.0;      // no WD Teff,
            setMags (mag, cmpnt, mass, pCluster, pStar, evoModels, filters, ltau, globalMags);

            deriveCombinedMags (mag, clusterAv, flux, pCluster, pStar, evoModels, filters);
            /* now have magnitudes, want posterior probability */
            tmpLogPost = logPost1Star (pStar, pCluster, evoModels, filterPriorMin, filterPriorMax);
            tmpLogPost += log (dMass);
            tmpLogPost += log ((isochrone.mass[i + 1] - isochrone.mass[i]) / mass[0]);
            tmpPost = exp (tmpLogPost);

            (*post) += tmpPost;

        }
    }
}
