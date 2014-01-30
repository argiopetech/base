#include <array>
#include <utility>
#include <vector>

#include <cmath>

#include "Cluster.hpp"
#include "Star.hpp"

#include "evolve.hpp"
#include "WdCoolingModel.hpp"
#include "WdAtmosphereModel.hpp"
#include "structures.hpp"
#include "Model.hpp"
#include "ifmr.hpp"

using std::array;
using std::vector;

double wdEvol (const Cluster &pCluster, const Model &evoModels, const vector<int> &filters, array<double, FILTS> &globalMags, Star &pStar, int cmpnt)
{
    std::pair<double, double> teffRadiusPair;

    double thisWDMass = 0.0, thisPrecLogAge = 0.0, thisLogTeff, thisWDLogRadius = 0.0;
    double mass;

    if (cmpnt == 0)
        mass = pStar.getMass1();
    else
        mass = pStar.getMass2();

    thisPrecLogAge = evoModels.mainSequenceEvol->wdPrecLogAge(pCluster.feh, mass);

    thisWDMass = intlFinalMassReln (pCluster, evoModels, mass);

    //get temperature from WD cooling models (returns 0.0 if there is an error(or does it??))
    teffRadiusPair = evoModels.WDcooling->wdMassToTeffAndRadius (pCluster.age, pCluster.carbonicity, thisPrecLogAge, thisWDMass);

    thisLogTeff = teffRadiusPair.first;
    thisWDLogRadius = teffRadiusPair.second;

    //*******this now gets trapped for in wdMassToTeffAndRadius so it should be unnecessary here (???)
    if (thisPrecLogAge >= pCluster.age)
    {                           // mcmc.c can cause this by adjusting masses and ages
        for (auto f : filters)
            globalMags[f] = -4.; // place at tip of RGB
    }
    else
    {
        //Calculate logg
        globalMags = evoModels.WDAtmosphere->teffToMags (thisLogTeff, thisWDMass, pStar.wdType[cmpnt]);
    }

    pStar.massNow[cmpnt] = thisWDMass;
    pStar.wdLogTeff[cmpnt] = thisLogTeff;
    pStar.status[cmpnt] = WD;

    return thisPrecLogAge;
}
