#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "evolve.h"
#include "linInterp.h"
#include "wdCooling.h"
#include "binSearch.h"
#include <mpi.h>

static int    nIso;
static double *wdMasses;
static struct wdCoolingCurve *wdCurves;
static int coolingModel;

void loadWDCool(char *path, int modelSet)
{
    int massCurves = 0, carbonCurves = 0, entries = 0;
    FILE *pCoolingModels;
    char tempFile[100]="\0", line[240];
    double tempAge, tempTeff, tempMass, tempRadius, lastMass=0.0, tempCarbon, lastCarbon = 0.0;
    void *tempAlloc; // temporary for allocation

    coolingModel = modelSet;

    // Allocate memory for the cooling curves dynamically so that
    // any number of curves can be read in (up to limits of memory)
    if ((wdCurves = (struct wdCoolingCurve *) calloc(1, sizeof(struct wdCoolingCurve))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");

    if ((wdCurves[0].wdCarbons = (double *) calloc(1, sizeof(double))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");

    if ((wdCurves[0].carbonCurve = (struct wdCarbonCurve *) calloc(1, sizeof(struct wdCarbonCurve))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");

    if ((wdMasses = (double *) calloc(1, sizeof(double))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");

    strcat(tempFile, path);

    if (modelSet == WOOD)
    {
        strcat(tempFile,"xb.comb");
    }
    else if (modelSet == MONTGOMERY)
    {
        strcat(tempFile,"wdtables");
    }
    else
    {
        printf("\nCooling models do not exist.  Exiting...\n");
        exit(1);
    }

    if ((pCoolingModels = fopen(tempFile,"r")) == NULL)
    {
        printf("\n\n file %s was not found - exiting\n",tempFile);
        exit(1);
    }

    fgets(line, 240, pCoolingModels); /* after header line, read in Wood model file */

    while (fgets(line, 240, pCoolingModels) != NULL)
    {
        if (modelSet == WOOD)
        {
            sscanf(line,"%*d %lf %*f %*f %*f %lf %lf %*f %*f %*f %lf", &tempAge, &tempRadius, &tempTeff, &tempMass);
            tempCarbon = 0.6; // Good extimate per 5 March, 2013 conversation with Dr. von Hippel
        }
        else
        {
            sscanf(line,"%*d %lf %*f %*f %*f %lf %lf %*f %*f %*f %lf %lf", &tempAge, &tempRadius, &tempTeff, &tempMass, &tempCarbon);
        }

        if ((massCurves == 0) && (carbonCurves == 0) && (entries == 0))
        {
            lastMass = tempMass;
            lastCarbon = tempCarbon;
        }

        // If the mass for this entry isn't the same as the mass of
        // the last entry, it's a new cooling curve.  Re-allocate
        // memory for the wdCurves array and start storing in the next
        // entry. Additionally, make a new cooling curve every time
        // there is a new carbon. In theory it will already do this
        // due to going from max to min mass, but I check anyway.

        if ((tempMass != lastMass))
        {
            wdCurves[massCurves].carbonCurve[carbonCurves].length = entries;
            wdCurves[massCurves].length = carbonCurves + 1;

            massCurves += 1;
            carbonCurves = 0;
            entries = 0;

            if ((tempAlloc = (void *) realloc(wdCurves, (massCurves+1) * sizeof(struct wdCoolingCurve))) == NULL)
                perror("wdCurves memory allocation error \n");
            else
                wdCurves = (struct wdCoolingCurve *) tempAlloc;

            if ((wdCurves[massCurves].wdCarbons = (double *) calloc(1, sizeof(double))) == NULL)
                perror("MEMORY ALLOCATION ERROR \n");

            if ((wdCurves[massCurves].carbonCurve = (struct wdCarbonCurve *) calloc(1, sizeof(struct wdCarbonCurve))) == NULL)
                perror("MEMORY ALLOCATION ERROR \n");

            if ((tempAlloc = (void *) realloc(wdMasses, (massCurves+1) * sizeof(double))) == NULL)
                perror("wdMasses memory allocation error \n");
            else
                wdMasses = (double *) tempAlloc;

        }
        else if (tempCarbon != lastCarbon)
        {
            wdCurves[massCurves].carbonCurve[carbonCurves].length = entries;
            carbonCurves += 1;
            entries = 0;

            if ((tempAlloc = (void *) realloc(wdCurves[massCurves].carbonCurve, (carbonCurves+1) * sizeof(struct wdCarbonCurve))) == NULL)
                perror("wdMasses memory allocation error \n");
            else
                wdCurves[massCurves].carbonCurve = (struct wdCarbonCurve *) tempAlloc;

            if ((tempAlloc = (void *) realloc(wdCurves[massCurves].wdCarbons, (carbonCurves+1) * sizeof(double))) == NULL)
                perror("wdMasses memory allocation error \n");
            else
                wdCurves[massCurves].wdCarbons = (double *) tempAlloc;
        }

        wdCurves[massCurves].carbonCurve[carbonCurves].logTeff[entries] = tempTeff;
        wdCurves[massCurves].mass = wdMasses[massCurves] = tempMass;
        wdCurves[massCurves].carbonCurve[carbonCurves].x_carbon = wdCurves[massCurves].wdCarbons[carbonCurves] = tempCarbon;
        wdCurves[massCurves].carbonCurve[carbonCurves].logRadius[entries] = tempRadius;
        wdCurves[massCurves].carbonCurve[carbonCurves].logAge[entries] = log10(tempAge);
        lastMass = tempMass;
        lastCarbon = tempCarbon;
        entries++;
    }

    wdCurves[massCurves].length = carbonCurves + 1;
    wdCurves[massCurves].carbonCurve[carbonCurves].length = entries;

    nIso = massCurves + 1;

    fclose(pCoolingModels);

    /* int i, j; */
    /* for (i = 0; i < nIso; i++) */
    /* { */
    /*     printf("%d\n", wdCurves[i].length); */
    /*     for (j = 0; j < wdCurves[i].length; j++) */
    /*     { */
    /*         printf("%g -- %g, %d\n", wdCurves[i].mass, wdCurves[i].carbonCurve[j].x_carbon, wdCurves[i].carbonCurve[j].length); */
    /*     } */
    /* } */

    /* exit(0); */
}


double wdMassToTeffAndRadius(double logAge, double x_carbon, double wdPrecLogAge, double wdMass, double *thisWDLogRadius)
{
    if (coolingModel == WOOD)
    {
        return wdMassToTeffAndRadius_wood(logAge, wdPrecLogAge, wdMass, thisWDLogRadius);
    }
    else
    {
        return wdMassToTeffAndRadius_montgomery(logAge, x_carbon, wdPrecLogAge, wdMass, thisWDLogRadius);
    }

}

double wdMassToTeffAndRadius_montgomery(double logAge, double x_carbon, double wdPrecLogAge, double wdMass, double *thisWDLogRadius)
{
    double wdCoolLogAge = 0.0;
    int massIndex = -1, m, c, ageIndex,carbonIndex;
    double ageTeff[2] = {0,0};
    double ageRadius[2] = {0,0};
    double newTeff[2] = {0, 0};
    double newRadius[2] = {0, 0};

    if (logAge > wdPrecLogAge)
    {                 // age limit check: otherwise wdCoolLogAge
        wdCoolLogAge = log10(pow(10.0, logAge) - pow(10.0, wdPrecLogAge));
    }
    else
    {                                // mcmc.c can cause this by adjusting masses and ages
        (*thisWDLogRadius) = 0.0;
        return 0.0;                  // no need to calculate anything, return to evolve.c here
    }

    massIndex = binarySearch(wdMasses, nIso, wdMass);

    carbonIndex = binarySearch(wdCurves[massIndex].wdCarbons, wdCurves[massIndex].length, x_carbon);

    if (massIndex < 0)
    {
        printf("Error in binary search on mass (wdCooling.c)\n");
        exit(1);
    }

    //For each mass entry, interpolate in age
    for (m = massIndex; m <= massIndex + 1; m++)
    {
        for (c = carbonIndex; c <= carbonIndex + 1; c++)
        {
            ageIndex = binarySearch( wdCurves[m].carbonCurve[c].logAge, wdCurves[m].carbonCurve[c].length, wdCoolLogAge );

            if (ageIndex < 0)
            {
                printf("Error in binary search on age (wdCooling.c)\n");
                exit(1);
            }

            ageTeff[c - carbonIndex]   = linInterpExtrap( wdCurves[m].carbonCurve[c].logAge[ageIndex]   , wdCurves[m].carbonCurve[c].logAge[ageIndex+1]
                                                        , wdCurves[m].carbonCurve[c].logTeff[ageIndex]  , wdCurves[m].carbonCurve[c].logTeff[ageIndex+1]  , wdCoolLogAge);

            ageRadius[c - carbonIndex] = linInterpExtrap( wdCurves[m].carbonCurve[c].logAge[ageIndex]   , wdCurves[m].carbonCurve[c].logAge[ageIndex+1]
                                                        , wdCurves[m].carbonCurve[c].logRadius[ageIndex], wdCurves[m].carbonCurve[c].logRadius[ageIndex+1], wdCoolLogAge);
        }

        // Now interpolate in mass
        newTeff[m - massIndex] = linInterpExtrap(wdCurves[m].wdCarbons[carbonIndex], wdCurves[m].wdCarbons[carbonIndex + 1], ageTeff[0], ageTeff[1], x_carbon);

        newRadius[m - massIndex]  = linInterpExtrap(wdCurves[m].wdCarbons[carbonIndex], wdCurves[m].wdCarbons[carbonIndex + 1], ageRadius[0], ageRadius[1], x_carbon);
    }

    (*thisWDLogRadius) = linInterpExtrap(wdMasses[massIndex], wdMasses[massIndex + 1], newRadius[0], newRadius[1], wdMass);

    return linInterpExtrap(wdMasses[massIndex], wdMasses[massIndex + 1], newTeff[0], newTeff[1], wdMass);
}


double wdMassToTeffAndRadius_wood(double logAge, double wdPrecLogAge, double wdMass, double *thisWDLogRadius)
{
    double wdCoolLogAge = 0.0, newTeff = 0.0;
    int massIndex = -1, m, ageIndex;
    double ageTeff[2] = {0,0};
    double ageRadius[2] = {0,0};

    if (logAge > wdPrecLogAge)
    {                 // age limit check: otherwise wdCoolLogAge
        wdCoolLogAge = log10(pow(10.0, logAge) - pow(10.0, wdPrecLogAge));
    }
    else
    {                                // mcmc.c can cause this by adjusting masses and ages
        (*thisWDLogRadius) = 0.0;
        return 0.0;                  // no need to calculate anything, return to evolve.c here
    }

    massIndex = binarySearch(wdMasses, nIso, wdMass);

    if (massIndex < 0)
    {
        printf("Error in binary search on mass (wdCooling.c)\n");
        exit(1);
    }

    //For each mass entry, interpolate in age
    for (m = massIndex; m <= massIndex + 1; m++)
    {
        ageIndex = binarySearch( wdCurves[m].carbonCurve[0].logAge, wdCurves[m].carbonCurve[0].length, wdCoolLogAge );

        if (ageIndex < 0)
        {
            printf("Error in binary search on age (wdCooling.c)\n");
            exit(1);
        }

        ageTeff[m - massIndex] = linInterpExtrap(wdCurves[m].carbonCurve[0].logAge[ageIndex],  wdCurves[m].carbonCurve[0].logAge[ageIndex+1], wdCurves[m].carbonCurve[0].logTeff[ageIndex], wdCurves[m].carbonCurve[0].logTeff[ageIndex+1], wdCoolLogAge);

        ageRadius[m - massIndex] = linInterpExtrap(wdCurves[m].carbonCurve[0].logAge[ageIndex],  wdCurves[m].carbonCurve[0].logAge[ageIndex+1], wdCurves[m].carbonCurve[0].logRadius[ageIndex], wdCurves[m].carbonCurve[0].logRadius[ageIndex+1], wdCoolLogAge);

    }

    // Now interpolate in mass
    newTeff = linInterpExtrap(wdMasses[massIndex], wdMasses[massIndex + 1], ageTeff[0], ageTeff[1], wdMass);
    (*thisWDLogRadius) = linInterpExtrap(wdMasses[massIndex], wdMasses[massIndex + 1], ageRadius[0], ageRadius[1], wdMass);

    return newTeff;

}
