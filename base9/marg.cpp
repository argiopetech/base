#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "evolve.hpp"
#include "msRgbEvol.hpp"
#include "gBaraffeMag.hpp"
#include "densities.hpp"

const int MAX_ENTRIES = 370;

extern int verbose, useFilt[FILTS], aFilt, needMassNow;
extern double ltau[2];

// Used by sub-methods of msRgbEvol (gGirMag, gChabMag, etc...) and wdEvol (gBergMag)
extern double globalMags[FILTS];
extern double ageLimit[2];

extern struct globalIso isochrone;
static double clusterAbs[FILTS] = { 0 };

double wdEvol (struct cluster *pCluster, struct star *pStar, int cmpnt);
void setMags (double mag[][FILTS], int cmpnt, double *mass, struct cluster *pCluster, struct star *pStar);
void deriveCombinedMags (double mag[][FILTS], double clusterAv, double *flux, struct cluster *pCluster, struct star *pStar);
void calcPost (double *post, double dMass, double mag[][FILTS], double clusterAv, double *flux, double *mass, struct cluster *pCluster, struct star *pStar);

void calcAbsCoeffsForMarg (int filterSet);

/* evaluate on a grid of primary mass and mass ratio to approximate
   the integral */
double margEvolveWithBinary (struct cluster *pCluster, struct star *pStar)
{
    double mag[3][FILTS], mass[2], flux, clusterAv;
    double post = 0.0;

    //Don't recalculate AGB mass (and isochrone) if these parameters are the same as they
    //were last time through
    deriveAgbTipMass (pCluster);        // determine AGBt ZAMS mass, to find evol state

    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (pCluster->AGBt_zmass < EPS)
    {
        pStar->boundsFlag = 1;
        return -HUGE_VAL;
    }

    clusterAv = getParameter (pCluster, ABS);
    if (fabs (clusterAbs[0]) < EPS)
        calcAbsCoeffsForMarg (pCluster->evoModels.filterSet);

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

                calcPost (&post, dMass, mag, clusterAv, &flux, mass, pCluster, pStar);
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

void setMags (double mag[][FILTS], int cmpnt, double *mass, struct cluster *pCluster, struct star *pStar)
{
    int filt;

    if (mass[cmpnt] <= 0.0001)
    {                           // for non-existent secondary stars
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = 99.999;
        pStar->status[cmpnt] = DNE;
        pStar->massNow[cmpnt] = 0.0;
    }
    else if (mass[cmpnt] <= pCluster->AGBt_zmass)
    {                           // for main seq or giant star
        pStar->massNow[cmpnt] = msRgbEvol (pCluster, mass[cmpnt]);
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = globalMags[filt];
        pStar->status[cmpnt] = MSRG;    // keep track of evolutionary state
    }
    else if (mass[cmpnt] <= pCluster->M_wd_up)
    {                           // for white dwarf
        ltau[cmpnt] = wdEvol (pCluster, pStar, cmpnt);
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = globalMags[filt];
    }
    else if (mass[cmpnt] <= 100.)
    {                           // for neutron star or black hole remnant
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = 99.999;
        pStar->status[cmpnt] = NSBH;
    }
    else
    {
        if (verbose)
            printf (" This condition should not happen, %.2f greater than 100 Mo\n", mass[cmpnt]);
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                mag[cmpnt][filt] = 99.999;
        pStar->status[cmpnt] = DNE;
    }
}

void deriveCombinedMags (double mag[][FILTS], double clusterAv, double *flux, struct cluster *pCluster, struct star *pStar)
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
            mag[2][filt] += getParameter (pCluster, MOD);
            mag[2][filt] += (clusterAbs[filt] - 1.0) * clusterAv;       // add A_[u-k] (standard defn of modulus already includes Av)
            pStar->photometry[i++] = mag[2][filt];
            //i++;
        }
    }
}


void calcPost (double *post, double dMass, double mag[][FILTS], double clusterAv, double *flux, double *mass, struct cluster *pCluster, struct star *pStar)
{
    setMass1 (pStar, pCluster, mass[0]);

    int cmpnt = 0, filt;

    setMags (mag, cmpnt, mass, pCluster, pStar);

    double tmpLogPost, tmpPost;

    /* first try 0.0 massRatio */
    cmpnt = 1;
    pStar->massRatio = 0.0;
    for (filt = 0; filt < FILTS; filt++)
        if (useFilt[filt])
            globalMags[filt] = 99.999;
    pStar->massNow[cmpnt] = 0.0;
    ltau[cmpnt] = 0.0;          // may not be a WD, so no precursor age,
    pStar->wdLogTeff[cmpnt] = 0.0;      // no WD Teff,
    setMags (mag, cmpnt, mass, pCluster, pStar);

    deriveCombinedMags (mag, clusterAv, flux, pCluster, pStar);
    tmpLogPost = logPost1Star (pStar, pCluster);
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
            diffLow = pow (10.0, ((pStar->obsPhot[obsFilt] - nSD * sqrt (pStar->variance[obsFilt])) / -2.5)) - pow (10.0, (mag[0][filt] / -2.5));
            diffUp = pow (10.0, ((pStar->obsPhot[obsFilt] + nSD * sqrt (pStar->variance[obsFilt])) / -2.5)) - pow (10.0, (mag[0][filt] / -2.5));
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
            pStar->massRatio = mass[0] / isochrone.mass[i];
            for (filt = 0; filt < FILTS; filt++)
                if (useFilt[filt])
                    globalMags[filt] = 99.999;
            pStar->massNow[cmpnt] = 0.0;
            ltau[cmpnt] = 0.0;  // may not be a WD, so no precursor age,
            pStar->wdLogTeff[cmpnt] = 0.0;      // no WD Teff,
            setMags (mag, cmpnt, mass, pCluster, pStar);

            deriveCombinedMags (mag, clusterAv, flux, pCluster, pStar);
            /* now have magnitudes, want posterior probability */
            tmpLogPost = logPost1Star (pStar, pCluster);
            tmpLogPost += log (dMass);
            tmpLogPost += log ((isochrone.mass[i + 1] - isochrone.mass[i]) / mass[0]);
            tmpPost = exp (tmpLogPost);

            (*post) += tmpPost;

        }
    }
}


void calcAbsCoeffsForMarg (int filterSet)
{
    if (filterSet == UBVRIJHK)
    {
        clusterAbs[0] = 1.569;  // Cardelli, Clayton, Mathis 1989, table 3
        clusterAbs[1] = 1.337;  // yields A_u -> A_k = f(A_v), for standard filters
        clusterAbs[2] = 1.0;
        clusterAbs[3] = 0.751;
        clusterAbs[4] = 0.479;
        clusterAbs[5] = 0.282;
        clusterAbs[6] = 0.190;
        clusterAbs[7] = 0.114;
    }
    else if (filterSet == SDSS)
    {
        clusterAbs[0] = 5.155 / 3.1;    // Stoughton et al. (2002, AJ, 123, 485)
        clusterAbs[1] = 3.793 / 3.1;    // Table 22, which gives Afilter/E(B-V )
        clusterAbs[2] = 2.751 / 3.1;    // We use Afilter/Av, so all are divided
        clusterAbs[3] = 2.086 / 3.1;    // by Rv = Av/E(B-V) = 3.1, consistent
        clusterAbs[4] = 1.479 / 3.1;    // Cardelli et al. (1989).
        clusterAbs[5] = 0.282;  // JHK come from Cardelli (see above)
        clusterAbs[6] = 0.190;
        clusterAbs[7] = 0.114;
    }
    else if (filterSet == ACS)
    {                           // from Table 14, Sirianni et al. (2005, PASP, 117, 1049)
        clusterAbs[0] = 4.081 / 3.1;    // they used R=3.1; also derived via Cardelli et al. values
        clusterAbs[1] = 3.634 / 3.1;    // Ext. Ratios A(P)/E(B-V) in ACS/WFC Filters for diff. SEDs
        clusterAbs[2] = 3.042 / 3.1;    // SED F435W F475W F550M F555W F606W F625W F775W F814W
        clusterAbs[3] = 3.177 / 3.1;    // O5  4.192 3.773 3.052 3.233 2.936 2.673 2.005 1.864
        clusterAbs[4] = 2.809 / 3.1;    // G2  4.081 3.634 3.042 3.177 2.809 2.637 1.982 1.825
        clusterAbs[5] = 2.637 / 3.1;    // M0  3.994 3.555 3.030 3.115 2.716 2.616 1.965 1.796
        clusterAbs[6] = 1.982 / 3.1;
        clusterAbs[7] = 1.825 / 3.1;    // using values for G2 star, usually ~2% of O5 or M0 value
    }
    else
    {
        printf ("filterSet %d not found.  Exiting. (evolve.c)\n", filterSet);
        exit (1);
    }
}
