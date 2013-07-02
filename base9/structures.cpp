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


double getMass1 (Star &pStar, Cluster &pCluster)
{
    double mass;

    mass = pStar.U + pStar.beta[AGE][0] * (pCluster.getAge() - pCluster.mean[AGE]) + pStar.beta[MOD][0] * (pCluster.getMod() - pCluster.mean[MOD]) + pStar.beta[FEH][0] * (pCluster.getFeH() - pCluster.mean[FEH]) + pStar.beta[YYY][0] * (pCluster.getY() - pCluster.mean[YYY]) + pStar.betaMassRatio[0] * pow (pStar.massRatio, pStar.betaMassRatio[1]);

    return mass;
}

double getMass2 (Star &pStar, Cluster &pCluster)
{
    double mass2;

    mass2 = getMass1 (pStar, pCluster) * pStar.massRatio;
    return mass2;
}

void setMass1 (Star &pStar, Cluster &pCluster, double newMass)
{
    pStar.U = newMass - (pStar.beta[AGE][0] * (pCluster.getAge() - pCluster.mean[AGE]) + pStar.beta[MOD][0] * (pCluster.getMod() - pCluster.mean[MOD]) + pStar.beta[FEH][0] * (pCluster.getFeH() - pCluster.mean[FEH]) + pStar.beta[YYY][0] * (pCluster.getY() - pCluster.mean[YYY]) + pStar.betaMassRatio[0] * pow (pStar.massRatio, pStar.betaMassRatio[1]));
}

void setMass2 (Star &pStar, Cluster &pCluster, double newMass)
{
    pStar.massRatio = newMass / getMass1 (pStar, pCluster);
}

//NOTE: nEntries and nFilts need to be set in the calling program before this will work properly.
void allocateGlobalIso (struct globalIso &newIso)
{
    if (newIso.nEntries == 0 || newIso.nFilts == 0)
    {
        cerr << "Cannot allocate memory for global isochrone with " << newIso.nEntries << " entries and/or " << newIso.nFilts << " filts" << endl;
        exit (1);
    }

    if (newIso.mass != NULL)
    {
        cerr << "globalIso memory already allocated" << endl;
        return;
    }

    newIso.mass = new double[newIso.nEntries]();
    newIso.massNow = new double[newIso.nEntries]();
    newIso.eep = new int[newIso.nEntries]();

    int j = 0;

    newIso.mag = new double*[newIso.nEntries]();

    for (j = 0; j < newIso.nEntries; j++)
    {
        newIso.mag[j] = new double[newIso.nFilts]();
    }
}


//NOTE: nEntries and nFilts need to be set in the calling program before this will work properly.
void freeGlobalIso (struct globalIso &newIso)
{
    delete[] (newIso.mass);
    delete[] (newIso.massNow);
    delete[] (newIso.eep);

    int j = 0;

    for (j = 0; j < newIso.nEntries; j++)
    {
        delete[] (newIso.mag[j]);
    }
    delete[] (newIso.mag);
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
