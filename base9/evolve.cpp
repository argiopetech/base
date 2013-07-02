#include <array>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "evolve.hpp"
#include "msRgbEvol.hpp"
#include "gBaraffeMag.hpp"
#include "wdEvol.hpp"
#include "FilterSet.hpp"

using std::array;
using std::vector;

extern int useFilt[FILTS], aFilt;

// Used by sub-methods of msRgbEvol (gGirMag, gChabMag, etc...) and wdEvol (gBergMag)
double globalMags[FILTS];
double ageLimit[2];

struct globalIso isochrone;

static array<double, FILTS> clusterAbs;

void calcAbsCoeffs (int filterSet);


/**************************************************************************************
last update: 25Aug10

This routine is the master subroutine of the simulation code and is organized so mcmc.c
can call this routine one pair of stars at a time for any of a range of hypothetical
stellar, cluster, or model properties.  This routine in turn calls other subroutines.
The parameters used by this routine are the cluster properties of age, metallicity, distance,
reddening, the model set (which indicates a combination of a stellar evolution model,
an initial-final mass relation, a WD cooling model, and a WD atmosphere model),  and the
ZAMS masses of the two stars.  This routine does not control binary fraction (parent
routines will take care of binaries by creating a primary and a secondary mass, the
latter of which can be zero), nor whether a WD is a DA or a DB (again, controlled by
parent routine).  Using all of these inputs, this routine updates the photometry of the
star structure in the whichever of the U through K filters are being used (stored in useFilt).

In order to facilitate the interface with the Yale-Yonsai Fortran code, this function has
been modified to accept an array of stars, an index, and the number of stars in the array.
If the index is negative, it will derive the photometry for all of the stars.  The parameter
numStars is the number of stars in the array.  This function has no bounds checking, so
NUMSTARS NEEDS TO BE INPUT CORRECTLY.  If the index is positive, it will use that element
of the stars array.  You can also feed it a pointer to a single star and an index of 0 to
get the photometry of a single star. -- SD
***************************************************************************************/
void evolve (Cluster &pCluster, Model const &evoModels, Star &star, array<double, 2> &ltau)
{
    int filt, i, cmpnt;
    double mag[3][FILTS], mass[2], flux, clusterAv;

    //Allocate memory to global isochrone(if it hasn't been allocated already)
    if (isochrone.mass == NULL)
    {
        isochrone.nEntries = 370;
        isochrone.nFilts = FILTS;
        allocateGlobalIso (isochrone);
    }

    //Don't recalculate AGB mass (and isochrone) if these parameters are the same as they
    //were last time through
    if (fabs (isochrone.FeH - pCluster.getFeH()) > EPS || fabs (isochrone.logAge - pCluster.getAge()) > EPS || fabs (isochrone.Y - pCluster.getY()) > EPS)
    {
        pCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(pCluster.getFeH(), pCluster.getY(), pCluster.getAge());    // determine AGBt ZAMS mass, to find evol state
    }

    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (pCluster.AGBt_zmass < EPS)
    {
        star.boundsFlag = 1;
        return;
    }

    clusterAv = pCluster.getAbs();
    if (fabs (clusterAbs[0]) < EPS)
        calcAbsCoeffs (evoModels.filterSet, clusterAbs);

    mass[0] = getMass1 (star, pCluster);
    mass[1] = getMass2 (star, pCluster);

    if (star.status[0] == BD)
    {
        getBaraffeMags (pCluster.getAge(), mass[0]);

        for (filt = 0; filt < 8; filt++)
            if (useFilt[filt])
                mag[2][filt] = 99.999;

        for (filt = 8; filt < FILTS; filt++)
        {
            if (useFilt[filt])
            {
                mag[2][filt] = globalMags[filt];
            }
        }
    }
    else
    {
        for (cmpnt = 0; cmpnt < 2; cmpnt++)
        {
            for (filt = 0; filt < FILTS; filt++)
                if (useFilt[filt])
                    globalMags[filt] = 99.999;

            star.massNow[cmpnt] = 0.0;
            ltau[cmpnt] = 0.0;      // may not be a WD, so no precursor age,
            star.wdLogTeff[cmpnt] = 0.0;        // no WD Teff,

            if (mass[cmpnt] <= 0.0001)
            {                       // for non-existent secondary stars
                for (filt = 0; filt < FILTS; filt++)
                    if (useFilt[filt])
                        mag[cmpnt][filt] = 99.999;
                star.status[cmpnt] = DNE;
                star.massNow[cmpnt] = 0.0;
            }
            else if (mass[cmpnt] <= pCluster.AGBt_zmass)
            {                       // for main seq or giant star
                star.massNow[cmpnt] = evoModels.mainSequenceEvol->msRgbEvol (mass[cmpnt]);
                for (filt = 0; filt < FILTS; filt++)
                    if (useFilt[filt])
                        mag[cmpnt][filt] = globalMags[filt];
                star.status[cmpnt] = MSRG;      // keep track of evolutionary state
            }
            else if (mass[cmpnt] <= pCluster.M_wd_up)
            {                       // for white dwarf
                ltau[cmpnt] = wdEvol (pCluster, evoModels, star, cmpnt);
                for (filt = 0; filt < FILTS; filt++)
                    if (useFilt[filt])
                        mag[cmpnt][filt] = globalMags[filt];
            }
            else if (mass[cmpnt] <= 100.)
            {                       // for neutron star or black hole remnant
                for (filt = 0; filt < FILTS; filt++)
                    if (useFilt[filt])
                        mag[cmpnt][filt] = 99.999;
                star.status[cmpnt] = NSBH;
            }
            else
            {
                //     log << (" This condition should not happen, %.2f greater than 100 Mo\n", mass[cmpnt]);
                for (filt = 0; filt < FILTS; filt++)
                    if (useFilt[filt])
                        mag[cmpnt][filt] = 99.999;
                star.status[cmpnt] = DNE;
            }
        }

        // can now derive combined mags
        if (mag[1][aFilt] < 99.)
        {                           // if there is a secondary star (aFilt set in parent program)
            for (filt = 0; filt < FILTS; filt++)
            {                       // (NOTE: useFilt shortcut may help once doing binaries)
                if (useFilt[filt])
                {
                    flux = pow (10.0, (mag[0][filt] / -2.5));       // add up the fluxes of the primary
                    flux += pow (10.0, (mag[1][filt] / -2.5));      // and the secondary
                    mag[2][filt] = -2.5 * log10 (flux);     // (these 3 lines take 5% of run time for N large)
                }                   // if primary mag = 99.999, then this works
            }
        }                           // to make the combined mag = secondary mag
        else
        {
            for (filt = 0; filt < FILTS; filt++)
                if (useFilt[filt])
                    mag[2][filt] = mag[0][filt];
        }
    }
    i = 0;
    for (filt = 0; filt < FILTS; filt++)
    {                               // can now add distance and absorption
        if (useFilt[filt])
        {
            mag[2][filt] += pCluster.getMod();
            mag[2][filt] += (clusterAbs[filt] - 1.0) * clusterAv;   // add A_[u-k] (standard defn of modulus already includes Av)
            star.photometry[i++] = mag[2][filt];
        }
    }
}
