#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "evolve.hpp"
#include "gBergMag.hpp"
#include "wdCooling.hpp"
#include "structures.hpp"
#include "Model.hpp"
#include "ifmr.hpp"

const double LOG_G_PLUS_LOG_M_SUN = 26.12302173752;

extern int useFilt[FILTS];
extern double globalMags[FILTS];

double wdEvol (const Cluster &pCluster, const Model &evoModels, Star &pStar, int cmpnt)
{

    int filt;
    double thisWDMass = 0.0, thisPrecLogAge = 0.0, thisLogTeff, thisWDLogRadius = 0.0, thisWDLogG = 0.0;
    double mass;

    if (cmpnt == 0)
        mass = pStar.getMass1(pCluster);
    else
        mass = pStar.getMass2(pCluster);

    thisPrecLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(pCluster.getFeH(), pCluster.getY(), mass);

    thisWDMass = intlFinalMassReln (pCluster, evoModels, mass);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    thisLogTeff = evoModels.WDcooling->wdMassToTeffAndRadius (pCluster.getAge(), pCluster.carbonicity, thisPrecLogAge, thisWDMass, thisWDLogRadius);

    //*******this now gets trapped for in wdMassToTeffAndRadius so it should be unnecessary here (???)
    if (thisPrecLogAge >= pCluster.getAge())
    {                           // mcmc.c can cause this by adjusting masses and ages
        for (filt = 0; filt < FILTS; filt++)
            if (useFilt[filt])
                globalMags[filt] = -4.; // place at tip of RGB
    }
    else
    {
        //Calculate logg
        thisWDLogG = LOG_G_PLUS_LOG_M_SUN + log (thisWDMass) - 2 * thisWDLogRadius;
        bergeronTeffToMags (thisLogTeff, thisWDLogG, pStar.wdType[cmpnt]);
    }

    pStar.massNow[cmpnt] = thisWDMass;
    pStar.wdLogTeff[cmpnt] = thisLogTeff;
    pStar.status[cmpnt] = WD;

    return thisPrecLogAge;

}
