#include <vector>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "evolve.hpp"
#include "structures.hpp"
#include "Cluster.hpp"

using std::vector;

void initModels (struct model *models)
{
    models->evoModel = 0;
    models->brownDwarfEvol = 0;
    models->mainSequenceEvol = 0;
    models->IFMR = 0;
    models->WDcooling = 0;
    models->WDatm = 0;
    models->filterSet = 0;
    models->numFilts = 0;
};


double getMass1 (Star *pStar, Cluster *pCluster)
{
    double mass;

    mass = pStar->U + pStar->beta[AGE][0] * (pCluster->getAge() - pCluster->mean[AGE]) + pStar->beta[MOD][0] * (pCluster->getMod() - pCluster->mean[MOD]) + pStar->beta[FEH][0] * (pCluster->getFeH() - pCluster->mean[FEH]) + pStar->beta[YYY][0] * (pCluster->getY() - pCluster->mean[YYY]) + pStar->betaMassRatio[0] * pow (pStar->massRatio, pStar->betaMassRatio[1]);

    return mass;
}

double getMass2 (Star *pStar, Cluster *pCluster)
{
    double mass2;

    mass2 = getMass1 (pStar, pCluster) * pStar->massRatio;
    return mass2;
}

void setMass1 (Star *pStar, Cluster *pCluster, double newMass)
{
    pStar->U = newMass - (pStar->beta[AGE][0] * (pCluster->getAge() - pCluster->mean[AGE]) + pStar->beta[MOD][0] * (pCluster->getMod() - pCluster->mean[MOD]) + pStar->beta[FEH][0] * (pCluster->getFeH() - pCluster->mean[FEH]) + pStar->beta[YYY][0] * (pCluster->getY() - pCluster->mean[YYY]) + pStar->betaMassRatio[0] * pow (pStar->massRatio, pStar->betaMassRatio[1]));
}

void setMass2 (Star *pStar, Cluster *pCluster, double newMass)
{
    pStar->massRatio = newMass / getMass1 (pStar, pCluster);
}

char filterNames[FILTS][10];

void setFilterNames (int filterSet)
{
    int f;

    for (f = 0; f < FILTS; f++)
        strcpy (filterNames[f], "\0");
    if (filterSet == UBVRIJHK)
    {
        strcat (filterNames[0], "U");
        strcat (filterNames[1], "B");
        strcat (filterNames[2], "V");
        strcat (filterNames[3], "R");
        strcat (filterNames[4], "I");
        strcat (filterNames[5], "J");
        strcat (filterNames[6], "H");
        strcat (filterNames[7], "K");
    }
    else if (filterSet == ACS)
    {
        strcat (filterNames[0], "F435W");
        strcat (filterNames[1], "F475W");
        strcat (filterNames[2], "F550M");
        strcat (filterNames[3], "F555W");
        strcat (filterNames[4], "F606W");
        strcat (filterNames[5], "F625W");
        strcat (filterNames[6], "F775W");
        strcat (filterNames[7], "F814W");
    }
    else if (filterSet == SDSS)
    {
        strcat (filterNames[0], "u");
        strcat (filterNames[1], "g");
        strcat (filterNames[2], "r");
        strcat (filterNames[3], "i");
        strcat (filterNames[4], "z");
        strcat (filterNames[5], "J");
        strcat (filterNames[6], "H");
        strcat (filterNames[7], "K");
    }
    else
    {
        printf ("\nfilterSet %d not available.  Exiting.\n", filterSet);
        exit (1);
    }

    strcat (filterNames[8], "IrB");
    strcat (filterNames[9], "IrR");
    strcat (filterNames[10], "SpB1");
    strcat (filterNames[11], "SpB2");
    strcat (filterNames[12], "SpB3");
    strcat (filterNames[13], "SpB4");
}

char *getFilterName (int index)
{
    return filterNames[index];
}

//NOTE: nEntries and nFilts need to be set in the calling program before this will work properly.
void allocateGlobalIso (struct globalIso *newIso)
{
    if (newIso->nEntries == 0 || newIso->nFilts == 0)
    {
        printf ("Cannot allocate memory for global isochrone with %d entries and/or %d filts\n", newIso->nEntries, newIso->nFilts);
        exit (1);
    }

    if (newIso->mass != NULL)
    {
        printf ("globalIso memory already allocated\n");
        return;
    }

    newIso->mass = new double[newIso->nEntries]();
    newIso->massNow = new double[newIso->nEntries]();
    newIso->eep = new int[newIso->nEntries]();

    int j = 0;

    newIso->mag = new double*[newIso->nEntries]();

    for (j = 0; j < newIso->nEntries; j++)
    {
        newIso->mag[j] = new double[newIso->nFilts]();
    }
}


//NOTE: nEntries and nFilts need to be set in the calling program before this will work properly.
void freeGlobalIso (struct globalIso *newIso)
{
    delete[] (newIso->mass);
    delete[] (newIso->massNow);
    delete[] (newIso->eep);

    int j = 0;

    for (j = 0; j < newIso->nEntries; j++)
    {
        delete[] (newIso->mag[j]);
    }
    delete[] (newIso->mag);
}


// Swaps two mass entries in a global isochrone (n and n-1)
void swapGlobalEntries (struct globalIso *thisIso, int n, int useFilt[FILTS])
{
    int filt, tempEep;
    double tempMass;
    vector<double> tempMag(thisIso->nFilts);

    tempMass = thisIso->mass[n];
    tempEep = thisIso->eep[n];
    for (filt = 0; filt < thisIso->nFilts; filt++)
        if (useFilt[filt])
            tempMag[filt] = thisIso->mag[n][filt];

    thisIso->mass[n] = thisIso->mass[n - 1];
    thisIso->eep[n] = thisIso->eep[n - 1];
    for (filt = 0; filt < thisIso->nFilts; filt++)
        if (useFilt[filt])
            thisIso->mag[n][filt] = thisIso->mag[n - 1][filt];

    thisIso->mass[n - 1] = tempMass;
    thisIso->eep[n - 1] = tempEep;
    for (filt = 0; filt < thisIso->nFilts; filt++)
        if (useFilt[filt])
            thisIso->mag[n - 1][filt] = tempMag[filt];
}
