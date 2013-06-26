#include <string>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>
#include <mpi.h>

#include "evolve.hpp"
#include "linInterp.hpp"
#include "wdCooling.hpp"
#include "binSearch.hpp"

using std::string;

static int nIso;
static double *wdMasses;
static struct wdCoolingCurve *wdCurves;
static int coolingModel;

struct althausModel
{
    const char *filename;
    int hasHLum;
    double mass;
};

struct renedoModel
{
    const char *filename;
    double mass;
};

void loadWDCool (string path, int modelSet)
{
    static struct althausModel althaus[] = {
        {"T045_1E4.Z0", false, 0.45},
        {"T047_1E4.Z0", false, 0.47},
        {"T05_1E4.Z0", false, 0.50},
        {"T052_1E4.Z0", false, 0.52},
        {"T054_1E4.Z0", false, 0.54},
        {"T056_1E4.Z0", false, 0.56},
        {"T058_1E4.Z0", true, 0.58},
        {"T06_1E4.Z0", true, 0.60},
        {"T062_1E4.Z0", true, 0.62},
        {"T064_1E4.Z0", true, 0.64},
        {"T066_1E4.Z0", true, 0.66},
        {"T068_1E4.Z0", true, 0.68},
        {"T07_1E4.Z0", true, 0.70},
        {"T072_1E4.Z0", true, 0.72},
        {"T074_1E4.Z0", true, 0.74},
        {"T076_1E4.Z0", true, 0.76},
        {"T078_1E4.Z0", true, 0.78},
        {"T08_1E4.Z0", true, 0.80},
        {"T082_1E4.Z0", true, 0.82},
        {"T084_1E4.Z0", true, 0.84},
        {"T09_1E4.Z0", true, 0.90},
        {"T10_1E4.Z0", true, 1.00},
        {"T11_1E4.Z0", true, 1.10},
        {0, 0, 0}
    };

    static struct renedoModel renedo[] = {
        {"wd0524_z001.trk", 0.524},
        {"wd0570_z001.trk", 0.570},
        {"wd0593_z001.trk", 0.593},
        {"wd0609_z001.trk", 0.609},
        {"wd0632_z001.trk", 0.632},
        {"wd0659_z001.trk", 0.659},
        {"wd0705_z001.trk", 0.705},
        {"wd0767_z001.trk", 0.767},
        {"wd0837_z001.trk", 0.837},
        {"wd0877_z001.trk", 0.877},
        {"wd0934_z001.trk", 0.934},
        {0, 0}
    };

    int massCurves = 0, carbonCurves = 0, entries = 0;
    FILE *pCoolingModels;
    char tempFile[100] = "\0", line[240];
    double tempAge, tempTeff, tempMass, tempRadius, lastMass = 0.0, tempCarbon, lastCarbon = 0.0;
    void *tempAlloc;            // temporary for allocation

    coolingModel = modelSet;

    // Allocate memory for the cooling curves dynamically so that
    // any number of curves can be read in (up to limits of memory)
    wdCurves = new struct wdCoolingCurve[1]();
    wdCurves[0].wdCarbons = new double[1]();
    wdCurves[0].carbonCurve = new struct wdCarbonCurve[1]();
    wdMasses = new double[1]();

    strcat (tempFile, path.c_str());

    if (modelSet == WOOD)
    {
        strcat (tempFile, "xb.comb");
    }
    else if (modelSet == MONTGOMERY)
    {
        strcat (tempFile, "wdtables");
    }
    else if ((modelSet != ALTHAUS) && (modelSet != RENEDO))
    {
        printf ("\nCooling models do not exist.  Exiting...\n");
        exit (1);
    }

    if ((modelSet == MONTGOMERY) || (modelSet == WOOD))
    {

        if ((pCoolingModels = fopen (tempFile, "r")) == NULL)
        {
            printf ("\n\n file %s was not found - exiting\n", tempFile);
            exit (1);
        }

        fgets (line, 240, pCoolingModels);      /* after header line, read in Wood model file */

        while (fgets (line, 240, pCoolingModels) != NULL)
        {
            if (modelSet == WOOD)
            {
                sscanf (line, "%*d %lf %*f %*f %*f %lf %lf %*f %*f %*f %lf", &tempAge, &tempRadius, &tempTeff, &tempMass);
                tempCarbon = 0.6;       // Good extimate per 5 March, 2013 conversation with Dr. von Hippel
            }
            else
            {
                sscanf (line, "%*d %lf %*f %*f %*f %lf %lf %*f %*f %*f %lf %lf", &tempAge, &tempRadius, &tempTeff, &tempMass, &tempCarbon);
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

                if ((tempAlloc = (void *) realloc (wdCurves, (massCurves + 1) * sizeof (struct wdCoolingCurve))) == NULL)
                    perror ("wdCurves memory allocation error \n");
                else
                    wdCurves = (struct wdCoolingCurve *) tempAlloc;

                wdCurves[massCurves].wdCarbons = new double[1]();
                wdCurves[massCurves].carbonCurve = new struct wdCarbonCurve[1]();

                if ((tempAlloc = (void *) realloc (wdMasses, (massCurves + 1) * sizeof (double))) == NULL)
                    perror ("wdMasses memory allocation error \n");
                else
                    wdMasses = (double *) tempAlloc;

            }
            else if (tempCarbon != lastCarbon)
            {
                wdCurves[massCurves].carbonCurve[carbonCurves].length = entries;
                carbonCurves += 1;
                entries = 0;

                if ((tempAlloc = (void *) realloc (wdCurves[massCurves].carbonCurve, (carbonCurves + 1) * sizeof (struct wdCarbonCurve))) == NULL)
                    perror ("wdMasses memory allocation error \n");
                else
                    wdCurves[massCurves].carbonCurve = (struct wdCarbonCurve *) tempAlloc;

                if ((tempAlloc = (void *) realloc (wdCurves[massCurves].wdCarbons, (carbonCurves + 1) * sizeof (double))) == NULL)
                    perror ("wdMasses memory allocation error \n");
                else
                    wdCurves[massCurves].wdCarbons = (double *) tempAlloc;
            }

            wdCurves[massCurves].carbonCurve[carbonCurves].logTeff[entries] = tempTeff;
            wdCurves[massCurves].mass = wdMasses[massCurves] = tempMass;
            wdCurves[massCurves].carbonCurve[carbonCurves].x_carbon = wdCurves[massCurves].wdCarbons[carbonCurves] = tempCarbon;
            wdCurves[massCurves].carbonCurve[carbonCurves].logRadius[entries] = tempRadius;
            wdCurves[massCurves].carbonCurve[carbonCurves].logAge[entries] = log10 (tempAge);
            lastMass = tempMass;
            lastCarbon = tempCarbon;
            entries++;
        }
    }
    else if (modelSet == ALTHAUS)
    {
        int i = 0;

        tempCarbon = 0.6;               // Good extimate per 5 March, 2013 conversation with Dr. von Hippel
        massCurves = -1;

        while (althaus[i].filename != 0)        // Keep going till we hit the last record
        {
            strcpy (tempFile, path.c_str());
            strcat (tempFile, "althaus/");
            strcat (tempFile, althaus[i].filename);

            if ((pCoolingModels = fopen (tempFile, "r")) == NULL)
            {
                printf ("\n\n file %s was not found - exiting\n", tempFile);
                exit (1);
            }

            massCurves += 1;
            carbonCurves = 0;
            entries = 0;

            if ((tempAlloc = (void *) realloc (wdCurves, (massCurves + 1) * sizeof (struct wdCoolingCurve))) == NULL)
                perror ("wdCurves memory allocation error \n");
            else
                wdCurves = (struct wdCoolingCurve *) tempAlloc;

            wdCurves[massCurves].wdCarbons = new double[1]();
            wdCurves[massCurves].carbonCurve = new struct wdCarbonCurve[1]();

            if ((tempAlloc = (void *) realloc (wdMasses, (massCurves + 1) * sizeof (double))) == NULL)
                perror ("wdMasses memory allocation error \n");
            else
                wdMasses = (double *) tempAlloc;

            fgets (line, 240, pCoolingModels);  // Read in two header lines
            fgets (line, 240, pCoolingModels);

            while (fgets (line, 240, pCoolingModels) != NULL)
            {
                if (althaus[i].hasHLum) // Has one extra throwaway field
                {
                    sscanf (line, "%*f %lf %*f %*f %*f %lf %lf %*f %*f %*f", &tempTeff, &tempAge, &tempRadius);
                }
                else
                {
                    sscanf (line, "%*f %lf %*f %*f %*f %lf %lf %*f %*f", &tempTeff, &tempAge, &tempRadius);
                }

                wdCurves[massCurves].carbonCurve[carbonCurves].logTeff[entries] = tempTeff;
                wdCurves[massCurves].mass = wdMasses[massCurves] = althaus[i].mass;
                wdCurves[massCurves].carbonCurve[carbonCurves].x_carbon = wdCurves[massCurves].wdCarbons[carbonCurves] = tempCarbon;
                wdCurves[massCurves].carbonCurve[carbonCurves].logRadius[entries] = tempRadius;
                wdCurves[massCurves].carbonCurve[carbonCurves].logAge[entries] = log10 (1e6) + tempAge;
                entries++;
                assert (entries < MAX_WD_MODEL);
            }

            wdCurves[massCurves].carbonCurve[carbonCurves].length = entries;
            wdCurves[massCurves].length = carbonCurves + 1;

//            fclose(pCoolingModels);
            pCoolingModels = 0;
            i += 1;
        }
    }
    else if (modelSet == RENEDO)
    {
        int i = 0;

        tempCarbon = 0.6;               // Good extimate per 5 March, 2013 conversation with Dr. von Hippel
        massCurves = -1;

        while (renedo[i].filename != 0) // Keep going till we hit the last record
        {
            strcpy (tempFile, path.c_str());
            strcat (tempFile, "renedo/");
            strcat (tempFile, renedo[i].filename);

            if ((pCoolingModels = fopen (tempFile, "r")) == NULL)
            {
                printf ("\n\n file %s was not found - exiting\n", tempFile);
                exit (1);
            }

            massCurves += 1;
            carbonCurves = 0;
            entries = 0;

            if ((tempAlloc = (void *) realloc (wdCurves, (massCurves + 1) * sizeof (struct wdCoolingCurve))) == NULL)
                perror ("wdCurves memory allocation error \n");
            else
                wdCurves = (struct wdCoolingCurve *) tempAlloc;


            wdCurves[massCurves].wdCarbons = new double[1]();
            wdCurves[massCurves].carbonCurve = new struct wdCarbonCurve[1]();

            if ((tempAlloc = (void *) realloc (wdMasses, (massCurves + 1) * sizeof (double))) == NULL)
                perror ("wdMasses memory allocation error \n");
            else
                wdMasses = (double *) tempAlloc;

            fgets (line, 240, pCoolingModels);  // Read in header line

            while (fgets (line, 240, pCoolingModels) != NULL)
            {
                sscanf (line, "%*f %lf %*f %*f %lf %*f %*f %*f %*f %*f %*f %*f %lf", &tempTeff, &tempAge, &tempRadius);

                if (tempAge >= 0)
                {
                    wdCurves[massCurves].carbonCurve[carbonCurves].logTeff[entries] = tempTeff;
                    wdCurves[massCurves].mass = wdMasses[massCurves] = renedo[i].mass;
                    wdCurves[massCurves].carbonCurve[carbonCurves].x_carbon = wdCurves[massCurves].wdCarbons[carbonCurves] = tempCarbon;
                    wdCurves[massCurves].carbonCurve[carbonCurves].logRadius[entries] = log10 (tempRadius * R_sun);
                    wdCurves[massCurves].carbonCurve[carbonCurves].logAge[entries] = log10 (1e6 * tempAge);
                    entries++;
                    assert (entries < MAX_WD_MODEL);
                }
            }

            if (tempAge >= 0)
            {
                wdCurves[massCurves].carbonCurve[carbonCurves].length = entries;
                wdCurves[massCurves].length = carbonCurves + 1;

                pCoolingModels = 0;
                i += 1;
            }
        }
    }

    wdCurves[massCurves].length = carbonCurves + 1;
    wdCurves[massCurves].carbonCurve[carbonCurves].length = entries;

    nIso = massCurves + 1;

//    fclose(pCoolingModels);

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


double wdMassToTeffAndRadius (double logAge, double x_carbon, double wdPrecLogAge, double wdMass, double *thisWDLogRadius)
{
    if ((coolingModel == WOOD) || (coolingModel == ALTHAUS) || (coolingModel == RENEDO))
    {
        return wdMassToTeffAndRadius_wood (logAge, wdPrecLogAge, wdMass, thisWDLogRadius);
    }
    else
    {
        return wdMassToTeffAndRadius_montgomery (logAge, x_carbon, wdPrecLogAge, wdMass, thisWDLogRadius);
    }

}

double wdMassToTeffAndRadius_montgomery (double logAge, double x_carbon, double wdPrecLogAge, double wdMass, double *thisWDLogRadius)
{
    double wdCoolLogAge = 0.0;
    int massIndex = -1, m, c, ageIndex, carbonIndex;
    double ageTeff[2] = { 0, 0 };
    double ageRadius[2] = { 0, 0 };
    double newTeff[2] = { 0, 0 };
    double newRadius[2] = { 0, 0 };

    if (logAge > wdPrecLogAge)
    {                           // age limit check: otherwise wdCoolLogAge
        wdCoolLogAge = log10 (pow (10.0, logAge) - pow (10.0, wdPrecLogAge));
    }
    else
    {                           // mcmc.c can cause this by adjusting masses and ages
        (*thisWDLogRadius) = 0.0;
        return 0.0;                     // no need to calculate anything, return to evolve.c here
    }

    massIndex = binarySearch (wdMasses, nIso, wdMass);

    carbonIndex = binarySearch (wdCurves[massIndex].wdCarbons, wdCurves[massIndex].length, x_carbon);

    if (massIndex < 0)
    {
        printf ("Error in binary search on mass (wdCooling.c)\n");
        exit (1);
    }

    //For each mass entry, interpolate in age
    for (m = massIndex; m <= massIndex + 1; m++)
    {
        for (c = carbonIndex; c <= carbonIndex + 1; c++)
        {
            ageIndex = binarySearch (wdCurves[m].carbonCurve[c].logAge, wdCurves[m].carbonCurve[c].length, wdCoolLogAge);

            if (ageIndex < 0)
            {
                printf ("Error in binary search on age (wdCooling.c)\n");
                exit (1);
            }

            ageTeff[c - carbonIndex] = linInterpExtrap (wdCurves[m].carbonCurve[c].logAge[ageIndex], wdCurves[m].carbonCurve[c].logAge[ageIndex + 1], wdCurves[m].carbonCurve[c].logTeff[ageIndex], wdCurves[m].carbonCurve[c].logTeff[ageIndex + 1], wdCoolLogAge);

            ageRadius[c - carbonIndex] = linInterpExtrap (wdCurves[m].carbonCurve[c].logAge[ageIndex], wdCurves[m].carbonCurve[c].logAge[ageIndex + 1], wdCurves[m].carbonCurve[c].logRadius[ageIndex], wdCurves[m].carbonCurve[c].logRadius[ageIndex + 1], wdCoolLogAge);
        }

        // Now interpolate in mass
        newTeff[m - massIndex] = linInterpExtrap (wdCurves[m].wdCarbons[carbonIndex], wdCurves[m].wdCarbons[carbonIndex + 1], ageTeff[0], ageTeff[1], x_carbon);

        newRadius[m - massIndex] = linInterpExtrap (wdCurves[m].wdCarbons[carbonIndex], wdCurves[m].wdCarbons[carbonIndex + 1], ageRadius[0], ageRadius[1], x_carbon);
    }

    (*thisWDLogRadius) = linInterpExtrap (wdMasses[massIndex], wdMasses[massIndex + 1], newRadius[0], newRadius[1], wdMass);

    return linInterpExtrap (wdMasses[massIndex], wdMasses[massIndex + 1], newTeff[0], newTeff[1], wdMass);
}


double wdMassToTeffAndRadius_wood (double logAge, double wdPrecLogAge, double wdMass, double *thisWDLogRadius)
{
    double wdCoolLogAge = 0.0, newTeff = 0.0;
    int massIndex = -1, m, ageIndex;
    double ageTeff[2] = { 0, 0 };
    double ageRadius[2] = { 0, 0 };

    if (logAge > wdPrecLogAge)
    {                           // age limit check: otherwise wdCoolLogAge
        wdCoolLogAge = log10 (pow (10.0, logAge) - pow (10.0, wdPrecLogAge));
    }
    else
    {                           // mcmc.c can cause this by adjusting masses and ages
        (*thisWDLogRadius) = 0.0;
        return 0.0;                     // no need to calculate anything, return to evolve.c here
    }

    massIndex = binarySearch (wdMasses, nIso, wdMass);

    if (massIndex < 0)
    {
        printf ("Error in binary search on mass (wdCooling.c)\n");
        exit (1);
    }

    //For each mass entry, interpolate in age
    for (m = massIndex; m <= massIndex + 1; m++)
    {
        ageIndex = binarySearch (wdCurves[m].carbonCurve[0].logAge, wdCurves[m].carbonCurve[0].length, wdCoolLogAge);

        if (ageIndex < 0)
        {
            printf ("Error in binary search on age (wdCooling.c)\n");
            exit (1);
        }

        ageTeff[m - massIndex] = linInterpExtrap (wdCurves[m].carbonCurve[0].logAge[ageIndex], wdCurves[m].carbonCurve[0].logAge[ageIndex + 1], wdCurves[m].carbonCurve[0].logTeff[ageIndex], wdCurves[m].carbonCurve[0].logTeff[ageIndex + 1], wdCoolLogAge);

        ageRadius[m - massIndex] = linInterpExtrap (wdCurves[m].carbonCurve[0].logAge[ageIndex], wdCurves[m].carbonCurve[0].logAge[ageIndex + 1], wdCurves[m].carbonCurve[0].logRadius[ageIndex], wdCurves[m].carbonCurve[0].logRadius[ageIndex + 1], wdCoolLogAge);

    }

    // Now interpolate in mass
    newTeff = linInterpExtrap (wdMasses[massIndex], wdMasses[massIndex + 1], ageTeff[0], ageTeff[1], wdMass);
    (*thisWDLogRadius) = linInterpExtrap (wdMasses[massIndex], wdMasses[massIndex + 1], ageRadius[0], ageRadius[1], wdMass);

    return newTeff;

}
