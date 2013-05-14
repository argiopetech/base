#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "evolve.h"
#include "msRgbEvol.h"
#include "gBaraffeMag.h"

extern int    verbose, useFilt[FILTS], aFilt, needMassNow;
extern double ltau[2];

// Used by sub-methods of msRgbEvol (gGirMag, gChabMag, etc...) and wdEvol (gBergMag)
double globalMags[FILTS];
double ageLimit[2];

struct globalIso isochrone;

static double clusterAbs[FILTS]={0};

double wdEvol(struct cluster *pCluster, struct star *pStar, int cmpnt);
void calcAbsCoeffs(int filterSet);

void evolve(struct cluster *pCluster, struct star *stars, int index)

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

{

  int    filt, i, cmpnt, j, numStars = 1;
  double mag[3][FILTS], mass[2], flux, clusterAv;

  //Allocate memory to global isochrone(if it hasn't been allocated already)
  if(isochrone.mass == NULL){
    isochrone.nEntries = 370;
    isochrone.nFilts = FILTS;
    allocateGlobalIso(&isochrone);
  }

  // A negative value for index means evolve all the stars in the *stars array,
  // a positive value means evolve that star only
  if(index < 0){
    index = 0;
    numStars = pCluster->nStars;
  }
  else numStars = 1;

  //Don't recalculate AGB mass (and isochrone) if these parameters are the same as they
  //were last time through
  if(fabs(isochrone.FeH - getParameter(pCluster,FEH)) > EPS ||
     fabs(isochrone.logAge - getParameter(pCluster,AGE)) > EPS ||
     fabs(isochrone.Y - getParameter(pCluster,YYY)) > EPS){
      deriveAgbTipMass(pCluster);                          // determine AGBt ZAMS mass, to find evol state
  }

  // AGBt_zmass never set because age and/or metallicity out of range of models.
  if(pCluster->AGBt_zmass < EPS){
    stars[0].boundsFlag = 1;
    return;
  }

  clusterAv   = getParameter(pCluster,ABS);
  if(fabs(clusterAbs[0]) < EPS) calcAbsCoeffs(pCluster->evoModels.filterSet);

  for(j = index; j < index + numStars; j++){
    mass[0]                   = getMass1(&stars[j],pCluster);
    mass[1]                   = getMass2(&stars[j],pCluster);

    if(stars[j].status[0] == BD){
      getBaraffeMags(getParameter(pCluster,AGE), mass[0]);
      for(filt=0;filt < 8;filt++) if(useFilt[filt]) mag[2][filt] = 99.999;
      for(filt=8;filt < FILTS;filt++){
        if(useFilt[filt]){
          mag[2][filt] = globalMags[filt];
        }
      }
    }
    else{
      for(cmpnt=0;cmpnt < 2;cmpnt++) {
        for(filt=0;filt < FILTS;filt++) if(useFilt[filt]) globalMags[filt] = 99.999;
        stars[j].massNow[cmpnt]   = 0.0;
        ltau[cmpnt]               = 0.0;                      // may not be a WD, so no precursor age,
        stars[j].wdLogTeff[cmpnt] = 0.0;                      // no WD Teff,

        if(mass[cmpnt] <= 0.0001) {                      // for non-existent secondary stars
          for(filt=0;filt < FILTS;filt++) if(useFilt[filt]) mag[cmpnt][filt] = 99.999;
          stars[j].status[cmpnt]  = DNE;
          stars[j].massNow[cmpnt] = 0.0;
        }
        else if(mass[cmpnt] <= pCluster->AGBt_zmass) {	// for main seq or giant star
          stars[j].massNow[cmpnt] = msRgbEvol(pCluster, mass[cmpnt]);
          for(filt=0;filt < FILTS;filt++) if(useFilt[filt]) mag[cmpnt][filt] = globalMags[filt];
          stars[j].status[cmpnt]  = MSRG;                  // keep track of evolutionary state
        }
        else if(mass[cmpnt] <= pCluster->M_wd_up) {	// for white dwarf
          ltau[cmpnt]   = wdEvol(pCluster, &(stars[j]), cmpnt);
          for(filt=0;filt < FILTS;filt++) if(useFilt[filt]) mag[cmpnt][filt] = globalMags[filt];
        }
        else if(mass[cmpnt] <= 100.) {                    // for neutron star or black hole remnant
          for(filt=0;filt < FILTS;filt++) if(useFilt[filt]) mag[cmpnt][filt] = 99.999;
          stars[j].status[cmpnt]  = NSBH;
        }
        else {
          if(verbose) printf(" This condition should not happen, %.2f greater than 100 Mo\n",mass[cmpnt]);
          for(filt=0;filt < FILTS;filt++) if(useFilt[filt]) mag[cmpnt][filt] = 99.999;
          stars[j].status[cmpnt]   = DNE;
        }
      }

      // can now derive combined mags
      if(mag[1][aFilt] < 99.) {                           // if there is a secondary star (aFilt set in parent program)
        for(filt=0;filt < FILTS;filt++) {                 // (NOTE: useFilt shortcut may help once doing binaries)
          if(useFilt[filt]){
            flux  = pow(10.0, (mag[0][filt] / -2.5));     // add up the fluxes of the primary
            flux += pow(10.0, (mag[1][filt] / -2.5));     // and the secondary
            mag[2][filt] = -2.5 * log10(flux);            // (these 3 lines take 5% of run time for N large)
          }                                               // if primary mag = 99.999, then this works
        }
      }                                                   // to make the combined mag = secondary mag
      else {
        for(filt=0;filt < FILTS;filt++) if(useFilt[filt]) mag[2][filt] = mag[0][filt];
      }
    }
    i=0;
    for(filt=0;filt < FILTS;filt++) {			// can now add distance and absorption
      if(useFilt[filt]) {
        mag[2][filt] += getParameter(pCluster,MOD);
        mag[2][filt] += (clusterAbs[filt] - 1.0)*clusterAv;	// add A_[u-k] (standard defn of modulus already includes Av)
        stars[j].photometry[i++] = mag[2][filt];
        //i++;
      }
    }
  }
}

void calcAbsCoeffs(int filterSet){
        if(filterSet == UBVRIJHK){
                clusterAbs[0] = 1.569;			// Cardelli, Clayton, Mathis 1989, table 3
                clusterAbs[1] = 1.337;			// yields A_u -> A_k = f(A_v), for standard filters
                clusterAbs[2] = 1.0;
                clusterAbs[3] = 0.751;
                clusterAbs[4] = 0.479;
                clusterAbs[5] = 0.282;
                clusterAbs[6] = 0.190;
                clusterAbs[7] = 0.114;
  }
        else if(filterSet == SDSS) {
                clusterAbs[0] = 5.155 / 3.1;            // Stoughton et al. (2002, AJ, 123, 485)
                clusterAbs[1] = 3.793 / 3.1;            // Table 22, which gives Afilter/E(B-V )
                clusterAbs[2] = 2.751 / 3.1;            // We use Afilter/Av, so all are divided
                clusterAbs[3] = 2.086 / 3.1;            // by Rv = Av/E(B-V) = 3.1, consistent
                clusterAbs[4] = 1.479 / 3.1;            // Cardelli et al. (1989).
                clusterAbs[5] = 0.282;                  // JHK come from Cardelli (see above)
                clusterAbs[6] = 0.190;
                clusterAbs[7] = 0.114;
        }
        else if(filterSet == ACS) {					// from Table 14, Sirianni et al. (2005, PASP, 117, 1049)
                clusterAbs[0] = 4.081 / 3.1;			// they used R=3.1; also derived via Cardelli et al. values
                clusterAbs[1] = 3.634 / 3.1;                    // Ext. Ratios A(P)/E(B-V) in ACS/WFC Filters for diff. SEDs
                clusterAbs[2] = 3.042 / 3.1;			// SED F435W F475W F550M F555W F606W F625W F775W F814W
                clusterAbs[3] = 3.177 / 3.1;			// O5  4.192 3.773 3.052 3.233 2.936 2.673 2.005 1.864
                clusterAbs[4] = 2.809 / 3.1;			// G2  4.081 3.634 3.042 3.177 2.809 2.637 1.982 1.825
                clusterAbs[5] = 2.637 / 3.1;			// M0  3.994 3.555 3.030 3.115 2.716 2.616 1.965 1.796
                clusterAbs[6] = 1.982 / 3.1;
                clusterAbs[7] = 1.825 / 3.1;			// using values for G2 star, usually ~2% of O5 or M0 value
        }
        else{
                printf("filterSet %d not found.  Exiting. (evolve.c)\n",filterSet);
                exit(1);
        }
}
