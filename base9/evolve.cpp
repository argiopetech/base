#include <array>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "evolve.hpp"
#include "marg.hpp"
#include "gBaraffeMag.hpp"
#include "wdEvol.hpp"
#include "FilterSet.hpp"
#include "WhiteDwarf.hpp"

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
void evolve (Cluster &pCluster, const Model &evoModels, Star &star, array<double, 2> &ltau)
{
    double mag[3][FILTS], mass[2], flux, clusterAv;

    //Don't recalculate AGB mass (and isochrone) if these parameters are the same as they
    //were last time through
    if (fabs (isochrone.FeH - pCluster.getFeH()) > EPS || fabs (isochrone.logAge - pCluster.getAge()) > EPS || fabs (isochrone.Y - pCluster.getY()) > EPS)
    {
        pCluster.AGBt_zmass = evoModels.mainSequenceEvol->deriveAgbTipMass(pCluster.getFeH(), pCluster.getY(), pCluster.getAge());    // determine AGBt ZAMS mass, to find evol state
    }

    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (pCluster.AGBt_zmass < EPS)
    {
        throw WDBoundsError("Bounds error in evolve.cpp");
    }

    clusterAv = pCluster.getAbs();

    if (fabs (clusterAbs[0]) < EPS)
        calcAbsCoeffs (evoModels.filterSet, clusterAbs);

    mass[0] = star.getMass1(pCluster);
    mass[1] = star.getMass2(pCluster);

    if (star.status[0] == BD)
    {
        getBaraffeMags (pCluster.getAge(), mass[0]);

        for (int filt = 0; filt < 8; filt++)
            if (useFilt[filt])
                mag[2][filt] = 99.999;

        for (int filt = 8; filt < FILTS; filt++)
        {
            if (useFilt[filt])
            {
                mag[2][filt] = globalMags[filt];
            }
        }
    }
    else
    {
        for (int cmpnt = 0; cmpnt < 2; cmpnt++)
        {
            setMags(mag, cmpnt, mass, pCluster, star, evoModels, ltau);
        }

        deriveCombinedMags(mag, clusterAv, &flux, pCluster, star);
    }
}
