#include <array>
#include <utility>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "evolve.hpp"
#include "gBergMag.hpp"
#include "WdCoolingModel.hpp"
#include "structures.hpp"
#include "Model.hpp"
#include "ifmr.hpp"

using std::array;
using std::vector;

const double LOG_G_PLUS_LOG_M_SUN = 26.12302173752;

double wdEvol (const Cluster &pCluster, const Model &evoModels, const vector<int> &filters, array<double, FILTS> &globalMags, Star &pStar, int cmpnt)
{
    std::pair<double, double> teffRadiusPair;

    double thisWDMass = 0.0, thisPrecLogAge = 0.0, thisLogTeff, thisWDLogRadius = 0.0, thisWDLogG = 0.0;
    double mass;

    if (cmpnt == 0)
        mass = pStar.getMass1();
    else
        mass = pStar.getMass2();

    thisPrecLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(pCluster.getFeH(), mass);

    thisWDMass = intlFinalMassReln (pCluster, evoModels, mass);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (pCluster.getAge(), pCluster.carbonicity, thisPrecLogAge, thisWDMass);

    thisLogTeff = teffRadiusPair.first;
    thisWDLogRadius = teffRadiusPair.second;

    //*******this now gets trapped for in wdMassToTeffAndRadius so it should be unnecessary here (???)
    if (thisPrecLogAge >= pCluster.getAge())
    {                           // mcmc.c can cause this by adjusting masses and ages
        for (auto f : filters)
            globalMags[f] = -4.; // place at tip of RGB
    }
    else
    {
        //Calculate logg
        thisWDLogG = LOG_G_PLUS_LOG_M_SUN + log (thisWDMass) - 2 * thisWDLogRadius;
        bergeronTeffToMags (filters, globalMags, thisLogTeff, thisWDLogG, pStar.wdType[cmpnt]);
    }

    pStar.massNow[cmpnt] = thisWDMass;
    pStar.wdLogTeff[cmpnt] = thisLogTeff;
    pStar.status[cmpnt] = WD;

    return thisPrecLogAge;

}
