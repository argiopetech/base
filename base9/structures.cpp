#include <iostream>
#include <vector>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "evolve.hpp"
#include "structures.hpp"
#include "Cluster.hpp"

using std::vector;
using std::cerr;
using std::endl;


double getMass1 (Star &pStar, const Cluster &pCluster)
{
    return pStar.U + pStar.beta[AGE][0] * (pCluster.getAge() - pCluster.mean[AGE])
                   + pStar.beta[MOD][0] * (pCluster.getMod() - pCluster.mean[MOD])
                   + pStar.beta[FEH][0] * (pCluster.getFeH() - pCluster.mean[FEH])
                   + pStar.beta[YYY][0] * (pCluster.getY() - pCluster.mean[YYY])
                   + pStar.betaMassRatio[0] * pow (pStar.massRatio, pStar.betaMassRatio[1]
                   );
}

double getMass2 (Star &pStar, const Cluster &pCluster)
{
    return getMass1 (pStar, pCluster) * pStar.massRatio;
}

void setMass1 (Star &pStar, const Cluster &pCluster, double newMass)
{
    pStar.U = newMass - ( pStar.beta[AGE][0] * (pCluster.getAge() - pCluster.mean[AGE]) 
                        + pStar.beta[MOD][0] * (pCluster.getMod() - pCluster.mean[MOD]) 
                        + pStar.beta[FEH][0] * (pCluster.getFeH() - pCluster.mean[FEH]) 
                        + pStar.beta[YYY][0] * (pCluster.getY()   - pCluster.mean[YYY])
                        + pStar.betaMassRatio[0] * pow (pStar.massRatio, pStar.betaMassRatio[1])
                        );
}

void setMass2 (Star &pStar, const Cluster &pCluster, double newMass)
{
    pStar.massRatio = newMass / getMass1 (pStar, pCluster);
}

// Swaps two mass entries in a global isochrone (n and n-1)
void swapGlobalEntries (struct globalIso &thisIso, int n, int useFilt[FILTS])
{
    int filt, tempEep;
    double tempMass;
    vector<double> tempMag(thisIso.nFilts);

    tempMass = thisIso.mass[n];
    tempEep = thisIso.eep[n];
    for (filt = 0; filt < thisIso.nFilts; filt++)
        if (useFilt[filt])
            tempMag[filt] = thisIso.mag[n][filt];

    thisIso.mass[n] = thisIso.mass[n - 1];
    thisIso.eep[n] = thisIso.eep[n - 1];
    for (filt = 0; filt < thisIso.nFilts; filt++)
        if (useFilt[filt])
            thisIso.mag[n][filt] = thisIso.mag[n - 1][filt];

    thisIso.mass[n - 1] = tempMass;
    thisIso.eep[n - 1] = tempEep;
    for (filt = 0; filt < thisIso.nFilts; filt++)
        if (useFilt[filt])
            thisIso.mag[n - 1][filt] = tempMag[filt];
}
