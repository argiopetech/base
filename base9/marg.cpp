#include <array>

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "evolve.hpp"
#include "msRgbEvol.hpp"
#include "gBaraffeMag.hpp"
#include "densities.hpp"
#include "Model.hpp"
#include "wdEvol.hpp"
#include "FilterSet.hpp"
#include "WhiteDwarf.hpp"

using std::array;

const int MAX_ENTRIES = 370;

extern int useFilt[FILTS], aFilt;

// Used by sub-methods of msRgbEvol (gGirMag, gChabMag, etc...) and wdEvol (gBergMag)
extern double globalMags[FILTS];
extern double ageLimit[2];

extern struct globalIso isochrone;
static array<double, FILTS> clusterAbs;

void setMags (double mag[][FILTS], int cmpnt, double *mass, Cluster &pCluster, Star &pStar, const Model&, array<double, 2>&);
void deriveCombinedMags (double mag[][FILTS], double clusterAv, double *flux, Cluster &pCluster, Star &pStar);
void calcPost (double *post, double dMass, double mag[][FILTS], double clusterAv, double *flux, double *mass, Cluster &pCluster, Star &pStar, const Model&, array<double, 2> &ltau);

/* evaluate on a grid of primary mass and mass ratio to approximate
   the integral */
double margEvolveWithBinary (Cluster &pCluster, Star &pStar, const Model &evoModels, array<double, 2> &ltau)
{
    double mag[3][FILTS], mass[2], flux, clusterAv;
    double post = 0.0;

    //Don't recalculate AGB mass (and isochrone) if these parameters are the same as they
    //were last time through
    pCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(pCluster.getFeH(), pCluster.getY(), pCluster.getAge());        // determine AGBt ZAMS mass, to find evol state

    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (pCluster.AGBt_zmass < EPS)
    {
        throw WDBoundsError("Bounds error in marg.cpp");
//        return -HUGE_VAL;
    }

    clusterAv = pCluster.getAbs();

    if (fabs (clusterAbs[0]) < EPS)
    {
        calcAbsCoeffs (evoModels.filterSet, clusterAbs);
    }

    int m;
    double dMass;

    double dIsoMass = 0.0;
    double k = 0.0;

    double isoIncrem = 80.0;    /* ok for YY models? */

    for (m = 0; m < isochrone.nEntries - 2; m++)
    {
        for (k = 0.0; k < isoIncrem; k += 1.0)
        {

            dIsoMass = isochrone.mass[m + 1] - isochrone.mass[m];

            /* why would dIsoMass ever be negative??? BUG in interpolation code??? */
            if (dIsoMass > 0.0)
            {
                dMass = dIsoMass / isoIncrem;
                mass[0] = isochrone.mass[m] + k * dMass;

                calcPost (&post, dMass, mag, clusterAv, &flux, mass, pCluster, pStar, evoModels, ltau);
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

void setMags (double mag[][FILTS], int cmpnt, double *mass, Cluster &pCluster, Star &pStar, const Model &evoModels, array<double, 2> &ltau)
{
    int filt;

    if (mass[cmpnt] <= 0.0001)
    {                           // for non-existent secondary stars
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = 99.999;
        pStar.status[cmpnt] = DNE;
        pStar.massNow[cmpnt] = 0.0;
    }
    else if (mass[cmpnt] <= pCluster.AGBt_zmass)
    {                           // for main seq or giant star
        pStar.massNow[cmpnt] = evoModels.mainSequenceEvol->msRgbEvol(mass[cmpnt]);
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = globalMags[filt];
        pStar.status[cmpnt] = MSRG;    // keep track of evolutionary state
    }
    else if (mass[cmpnt] <= pCluster.M_wd_up)
    {                           // for white dwarf
        ltau[cmpnt] = wdEvol (pCluster, evoModels, pStar, cmpnt);
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = globalMags[filt];
    }
    else if (mass[cmpnt] <= 100.)
    {                           // for neutron star or black hole remnant
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = 99.999;
        pStar.status[cmpnt] = NSBH;
    }
    else
    {
        //     log <<  (" This condition should not happen, %.2f greater than 100 Mo\n", mass[cmpnt]);
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = 99.999;
        pStar.status[cmpnt] = DNE;
    }
}

void deriveCombinedMags (double mag[][FILTS], double clusterAv, double *flux, Cluster &pCluster, Star &pStar)
{
    int filt;

    // can now derive combined mags
    if (mag[1][aFilt] < 99.)
    {                           // if there is a secondary star (aFilt set in parent program)
        for (filt = 0; filt < FILTS; filt++)
        {                               // (NOTE: useFilt shortcut may help once doing binaries)
            if (useFilt[filt])
            {
                (*flux) = pow (10.0, (mag[0][filt] / -2.5));    // add up the fluxes of the primary
                (*flux) += pow (10.0, (mag[1][filt] / -2.5));   // and the secondary
                mag[2][filt] = -2.5 * log10 (*flux);    // (these 3 lines take 5% of run time for N large)
            }                           // if primary mag = 99.999, then this works
        }
    }                           // to make the combined mag = secondary mag
    else
    {
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[2][filt] = mag[0][filt];
    }

    int i = 0;

    for (filt = 0; filt < FILTS; filt++)
    {                           // can now add distance and absorption
        if (useFilt[filt])
        {
            mag[2][filt] += pCluster.getMod();
            mag[2][filt] += (clusterAbs[filt] - 1.0) * clusterAv;       // add A_[u-k] (standard defn of modulus already includes Av)
            pStar.photometry[i++] = mag[2][filt];
            //i++;
        }
    }
}


void calcPost (double *post, double dMass, double mag[][FILTS], double clusterAv, double *flux, double *mass, Cluster &pCluster, Star &pStar, const Model &evoModels, array<double, 2> &ltau)
{
    setMass1 (pStar, pCluster, mass[0]);

    int cmpnt = 0, filt;

    setMags (mag, cmpnt, mass, pCluster, pStar, evoModels, ltau);

    double tmpLogPost, tmpPost;

    /* first try 0.0 massRatio */
    cmpnt = 1;
    pStar.massRatio = 0.0;
    for (filt = 0; filt < FILTS; filt++)
        if (useFilt[filt])
            globalMags[filt] = 99.999;
    pStar.massNow[cmpnt] = 0.0;
    ltau[cmpnt] = 0.0;          // may not be a WD, so no precursor age,
    pStar.wdLogTeff[cmpnt] = 0.0;      // no WD Teff,
    setMags (mag, cmpnt, mass, pCluster, pStar, evoModels, ltau);

    deriveCombinedMags (mag, clusterAv, flux, pCluster, pStar);
    tmpLogPost = logPost1Star (pStar, pCluster, evoModels);
    tmpLogPost += log (dMass);
    tmpLogPost += log (isochrone.mass[0] / mass[0]);    /* dMassRatio */
    tmpPost = exp (tmpLogPost);
    (*post) += tmpPost;


    /**** now see if any binary companions are OK ****/
    double magLower;
    double magUpper;
    double nSD = 4.0;           /* num of st dev from obs that will contribute to likelihood */
    int obsFilt = 0;
    int i, j;

    filt = 0;
    int isOverlap = 1;          /* do the allowable masses in each filter overlap? */
    double diffLow = 0.0;
    double diffUp = 0.0;
    int okMass[MAX_ENTRIES] = { 1 };
    while (filt < FILTS && isOverlap == 1)
    {
        if (useFilt[filt])
        {
            diffLow = pow (10.0, ((pStar.obsPhot[obsFilt] - nSD * sqrt (pStar.variance[obsFilt])) / -2.5)) - pow (10.0, (mag[0][filt] / -2.5));
            diffUp = pow (10.0, ((pStar.obsPhot[obsFilt] + nSD * sqrt (pStar.variance[obsFilt])) / -2.5)) - pow (10.0, (mag[0][filt] / -2.5));
            if (diffLow <= 0.0 || diffUp <= 0.0 || diffLow == diffUp)
            {
                isOverlap = 0;

                /**** necessary here??? *****/
                for (j = 0; j < FILTS; j++)
                    if (useFilt[j])
                        mag[2][j] = mag[0][j];
            }
            else
            {
                magLower = -2.5 * log10 (diffLow);
                magUpper = -2.5 * log10 (diffUp);

                for (i = 0; i < isochrone.nEntries - 1; i++)
                {
                    if (isochrone.mag[i][filt] >= magLower && isochrone.mag[i][filt] <= magUpper && isochrone.mass[i] <= mass[0])
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
        filt++;
    }

    for (i = 0; i < isochrone.nEntries - 2; i++)
    {
        if (okMass[i])
        {
            cmpnt = 1;
            pStar.massRatio = mass[0] / isochrone.mass[i];
            for (filt = 0; filt < FILTS; filt++)
                if (useFilt[filt])
                    globalMags[filt] = 99.999;
            pStar.massNow[cmpnt] = 0.0;
            ltau[cmpnt] = 0.0;  // may not be a WD, so no precursor age,
            pStar.wdLogTeff[cmpnt] = 0.0;      // no WD Teff,
            setMags (mag, cmpnt, mass, pCluster, pStar, evoModels, ltau);

            deriveCombinedMags (mag, clusterAv, flux, pCluster, pStar);
            /* now have magnitudes, want posterior probability */
            tmpLogPost = logPost1Star (pStar, pCluster, evoModels);
            tmpLogPost += log (dMass);
            tmpLogPost += log ((isochrone.mass[i + 1] - isochrone.mass[i]) / mass[0]);
            tmpPost = exp (tmpLogPost);

            (*post) += tmpPost;

        }
    }
}
