
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "evolve.hpp"
#include "gBergMag.hpp"
#include "wdCooling.hpp"
#include "structures.hpp"

#define LOG_G_PLUS_LOG_M_SUN 26.12302173752

extern int useFilt[FILTS];
extern double globalMags[FILTS];

double wdPrecLogAge (struct cluster *pCluster, double zamsMass);

// double intlFinalMassReln(double zamsMass, int IFMR);
double intlFinalMassReln (struct cluster *pCluster, double zamsMass);

double wdEvol (struct cluster *pCluster, struct star *pStar, int cmpnt)
/****************************************************************************************
last update: 02dec07

****************************************************************************************/
{

  int filt;
  double thisWDMass = 0.0, thisPrecLogAge = 0.0, thisLogTeff, thisWDLogRadius = 0.0, thisWDLogG = 0.0;
  double mass;

  if (cmpnt == 0)
    mass = getMass1 (pStar, pCluster);
  else
    mass = getMass2 (pStar, pCluster);

  thisPrecLogAge = wdPrecLogAge (pCluster, mass);
  // thisWDMass       = intlFinalMassReln(mass, pCluster->evoModels.IFMR);
  thisWDMass = intlFinalMassReln (pCluster, mass);

  //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
  // ***FIX ME***
  thisLogTeff = wdMassToTeffAndRadius (getParameter (pCluster, AGE), pCluster->carbonicity, thisPrecLogAge, thisWDMass, &thisWDLogRadius);

  //*******this now gets trapped for in wdMassToTeffAndRadius so it should be unnecessary here (???)
  if (thisPrecLogAge >= getParameter (pCluster, AGE))
  {				// mcmc.c can cause this by adjusting masses and ages
    for (filt = 0; filt < FILTS; filt++)
      if (useFilt[filt])
	globalMags[filt] = -4.;	// place at tip of RGB
  }
  else
  {
    //Calculate logg
    thisWDLogG = LOG_G_PLUS_LOG_M_SUN + log (thisWDMass) - 2 * thisWDLogRadius;
    bergeronTeffToMags (thisLogTeff, thisWDLogG, (*pStar).wdType[cmpnt]);
  }

  (*pStar).massNow[cmpnt] = thisWDMass;
  (*pStar).wdLogTeff[cmpnt] = thisLogTeff;
  (*pStar).status[cmpnt] = WD;

  return thisPrecLogAge;

}
