#include <string>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "evolve.hpp"
#include "ChabMsModel.hpp"
#include "binSearch.hpp"
#include "linInterp.hpp"

using std::string;

extern int verbose;
extern int useFilt[FILTS];
extern double globalMags[FILTS];
extern struct globalIso isochrone;

// defined in evolve.c, used to determine the normalization
// for the prior on age.
extern double ageLimit[2];

struct cIsochrone
{
    double age;                 //In Gyr
    double logAge;
    double FeH;
    double Y;
    int numEeps;
    int eeps[MAX_CHAB_ENTRIES];
    double mass[MAX_CHAB_ENTRIES];
    double mag[MAX_CHAB_ENTRIES][N_CHAB_FILTS];
    double AGBt;
};

static double cFeH[N_CHAB_Z], cY[N_CHAB_Z][N_CHAB_Y];
static double cLogAge[N_CHAB_Z][N_CHAB_Y][N_CHAB_AGES], cAge[N_CHAB_Z][N_CHAB_Y][N_CHAB_AGES];
static struct cIsochrone cIso[N_CHAB_Z][N_CHAB_Y][N_CHAB_AGES];

// Set in deriveChabAgbTip, used in getChaboyerMags
//static struct cIsochrone newIso;

static void calcCoeff (double a[], double b[], double x);

//static  void outputIso(struct cIsochrone *iso, FILE *wPtr);
static void initIso (struct cIsochrone *iso);
static void getFileName (string path, int z, int y);
static char tempFile[100];

void ChabMsModel::loadModel (string path, MsFilterSet filterSet)
{
    // ************************************************************************************
    // **********************    Chaboyer models   ****************************************
    // * BBN, dY/dZ=1.6 and Y_init=0.286, 0.306 consistent with initial solar composition *
    // ************************************************************************************


    int z, a, y, i;
    double tempMass = 0.0;
    char line[240];             //,tempFile[100];
    FILE *pChaboyer;

    if (filterSet != MsFilterSet::UBVRIJHK)
    {
        printf ("\nFilter set %d not available on Chaboyer helium models.  Exiting...\n", filterSet);
        exit (1);
    }

    for (z = 0; z < N_CHAB_Z; z++)
    {                           // foreach Chaboyer metallicity/isochrone file
        cFeH[z] = 0.0;

        cY[z][0] = 0.27;
        cY[z][1] = 0.30;
        cY[z][2] = 0.33;
        cY[z][3] = 0.36;
        cY[z][4] = 0.40;

        for (y = 0; y < N_CHAB_Y; y++)
        {                               // foreach Chaboyer metallicity/isochrone file

            for (a = 0; a < N_CHAB_AGES; a++)
            {                           // initialize age/boundary pointers
                cLogAge[z][y][a] = 0.0;
                initIso (&(cIso[z][y][a]));     // initialize array of model parameters
                cIso[z][y][a].Y = cY[z][y];
            }

            getFileName (path, z, y);

            //fscanf(pModelList,"%s",tempFile);                                 // work on one Chaboyer model at a time
            if ((pChaboyer = fopen (tempFile, "r")) == NULL)
            {                           // open file
                printf ("\n\n file %s was not found - exiting\n", tempFile);
                exit (1);
            }

            i = 0;
            a = -1;                     // a = [0,18], cLogAge contains the 19 ages
            while (fgets (line, 240, pChaboyer) != NULL)
            {                           // load each Z=# and Y=# Chaboyer model for all ages
                if (line[1] == 'M')
                {
                    fgets (line, 240, pChaboyer);
                    sscanf (line, "%*s %*f %*f %*f %*f %lf", &cFeH[z]);
                    if (fabs (cFeH[z]) < 0.02)
                        cFeH[z] = 0.0;
                }
                else if (line[1] == 'A')
                {
                    a++;
                    sscanf (line, "%*s %lf ", &(cIso[z][y][a].age));
                    cLogAge[z][y][a] = cIso[z][y][a].logAge = log10 (cIso[z][y][a].age * (1e9));
                    cAge[z][y][a] = cIso[z][y][a].age;
                    cIso[z][y][a].FeH = cFeH[z];
                    cIso[z][y][a].Y = cY[z][y];

                    if (a > 0)
                    {
                        cIso[z][y][a - 1].AGBt = tempMass;
                        cIso[z][y][a - 1].numEeps = i;
                    }
                    i = 0;
                }
                else if (line[0] != '#' && line[0] != '\n')
                {
                    sscanf (line, "%d %lf %*f %*f %*f %lf %lf %lf %lf %lf %lf %lf %lf", &cIso[z][y][a].eeps[i], &cIso[z][y][a].mass[i], &cIso[z][y][a].mag[i][0], &cIso[z][y][a].mag[i][1], &cIso[z][y][a].mag[i][2], &cIso[z][y][a].mag[i][3], &cIso[z][y][a].mag[i][4], &cIso[z][y][a].mag[i][5], &cIso[z][y][a].mag[i][6], &cIso[z][y][a].mag[i][7]);
                    tempMass = cIso[z][y][a].mass[i];
                    i++;
                }
            }
            // Put the last two entries into place
            cIso[z][y][a - 1].AGBt = tempMass;
            cIso[z][y][a - 1].numEeps = i;
        }
    }
    ageLimit[0] = cLogAge[0][0][0];
    ageLimit[1] = cLogAge[0][0][N_CHAB_AGES - 1];

    /*
      for(z=0 ; z < N_CHAB_Z ; z++) {                                       // foreach Chaboyer metallicity/isochrone file
      for(y=0 ; y < N_CHAB_Y ; y++) {                                     // foreach Chaboyer metallicity/isochrone file
      for(a=0 ; a < N_CHAB_AGES ; a++) {                               // initialize age/boundary pointers
      outputIso(&cIso[z][y][a], stdout);
      }
      }
      }
    */
}



static void getFileName (string path, int z, int y)
{
    char zString[][2] = { "0", "2", "4", "5" };
    char yString[][3] = { "27", "30", "33", "36", "39" };
    strcpy (tempFile, "\0");
    strcat (tempFile, path.c_str());
    strcat (tempFile, "cIso/feh0");
    strcat (tempFile, zString[z]);
    strcat (tempFile, "y\0");
    strcat (tempFile, yString[y]);
    strcat (tempFile, ".cmd\0");
}


// Interpolates between isochrones for two ages using linear interpolation
// Must run loadChaboyer() first for this to work.
double ChabMsModel::deriveAgbTipMass (double newFeH, double newY, double newLogAge)
{

    int iAge = -1, iY = -1, iFeH = -1, newimax = 500, newimin = 0, ioff[2][2][2], neweep;
    int z = 0, y = 0, a = 0, m = 0, filt = 0, n = 0;
    double newAge = pow (10, newLogAge) / 1e9;
    double b[2], c[2], d[2];

    //verbose = 1;
    //printf("******** Here  %f %f %f **********\n", newFeH, newY, newLogAge);

    if (newLogAge < cLogAge[0][0][0])
    {
        if (verbose)
            printf ("\n Requested age (%.3f) too young. (gChabMag.c)", cLogAge[0][0][0]);       //newLogAge);
        return 0.0;
    }
    if (newLogAge > cLogAge[N_CHAB_Z - 1][N_CHAB_Y - 1][N_CHAB_AGES - 1])
    {
        if (verbose)
            printf ("\n Requested age (%.3f) too old. (gChabMag.c)", newLogAge);
        return 0.0;
    }
    if (newFeH < cFeH[0])
    {
        if (verbose)
            printf ("\n Requested FeH (%.3f) too low. (gChabMag.c)", newFeH);
        return 0.0;
    }
    if (newFeH > cFeH[N_CHAB_Z - 1])
    {
        if (verbose)
            printf ("\n Requested FeH (%.3f) too high. (gChabMag.c)", newFeH);
        return 0.0;
    }
    if (newY < cY[0][0])
    {
        if (verbose)
            printf ("\n Requested Y (%.3f) too low. (gChabMag.c)", newY);
        return 0.0;
    }
    if (newY > cY[N_CHAB_Z - 1][N_CHAB_Y - 1])
    {
        if (verbose)
            printf ("\n Requested Y (%.3f) too high. (gChabMag.c)", newY);
        return 0.0;
    }

    // Find the values for each parameter that we will be interpolating
    // between and calculate the interpolation coefficients.
    iFeH = binarySearch (cFeH, N_CHAB_Z, newFeH);
    iY = binarySearch (cY[N_CHAB_Z - 1], N_CHAB_Y, newY);
    iAge = binarySearch (cAge[z][y], N_CHAB_AGES, newAge);

    calcCoeff (&cFeH[iFeH], d, newFeH);
    calcCoeff (&cY[z][iY], c, newY);
    calcCoeff (&cAge[z][y][iAge], b, newAge);

    // if(o) fprintf(wPtr,"%f %f %f %f %f %f\n",b[0],b[1],c[0],c[1],d[0],d[1]);

    // Find the minimum and maximum eep values and set a global
    // min and max that is the union of the min and max for each isochrone
    for (a = iAge; a < iAge + 2; a++)
    {
        for (y = iY; y < iY + 2; y++)
        {
            for (z = iFeH; z < iFeH + 2; z++)
            {
                if (cIso[z][y][a].eeps[0] > newimin)
                    newimin = cIso[z][y][a].eeps[0];
                if (cIso[z][y][a].eeps[cIso[z][y][a].numEeps - 1] < newimax)
                    newimax = cIso[z][y][a].eeps[cIso[z][y][a].numEeps - 1];
            }
        }
    }
    neweep = newimax - newimin + 1;     // = the # of eep points that will be in the new isochrone

    // For each isochrone, find the amount by which
    // the min eep # for that isochrone is
    // offset from the global minimum eep #
    for (a = 0; a < 2; a++)
    {
        for (y = 0; y < 2; y++)
        {
            for (z = 0; z < 2; z++)
            {
                ioff[z][y][a] = newimin - cIso[z + iFeH][y + iY][a + iAge].eeps[0];
            }
        }
    }

    // Now for each entry in the new isochrone
    // use the coefficients to calculate the mass
    // and each color at that eep
    for (m = 0; m < neweep; m++)
    {
        isochrone.mass[m] = 0.0;
        for (filt = 0; filt < N_CHAB_FILTS; filt++)
            if (useFilt[filt])
                isochrone.mag[m][filt] = 0.0;

        for (a = 0; a < 2; a++)
        {
            for (y = 0; y < 2; y++)
            {
                for (z = 0; z < 2; z++)
                {
                    isochrone.mass[m] += b[a] * c[y] * d[z] * cIso[iFeH + z][iY + y][iAge + a].mass[m + ioff[z][y][a]];
                    for (filt = 0; filt < N_CHAB_FILTS; filt++)
                    {
                        if (useFilt[filt])
                        {
                            isochrone.mag[m][filt] += b[a] * c[y] * d[z] * cIso[iFeH + z][iY + y][iAge + a].mag[m + ioff[z][y][a]][filt];
                        }
                    }
                }
            }
        }

        isochrone.eep[m] = cIso[iFeH][iY][iAge].eeps[m + ioff[0][0][0]];

        // Sometimes the interpolation process can leave the
        // mass entries out of order.  This swaps them so that
        // the later mass interpolation can function properly
        if (m > 0)
        {
            n = m;
            while (isochrone.mass[n] < isochrone.mass[n - 1] && n > 0)
            {
                swapGlobalEntries (&isochrone, n, useFilt);
                n--;
            }
        }
    }

    isochrone.nEntries = neweep;
    isochrone.age = newAge;
    isochrone.logAge = newLogAge;
    isochrone.FeH = newFeH;
    isochrone.Y = newY;
    isochrone.AgbTurnoffMass = isochrone.mass[neweep - 1];

    return isochrone.AgbTurnoffMass;

}


// Calculates magnitudes for a given mass.
// Must run loadChaboyer() and deriveChabAgbTip()
// to load and interpolate an isochrone before this subroutine will work
// Stores output values in globalMags[]
double ChabMsModel::msRgbEvol (double zamsMass)
{

    int m, filt;

    m = binarySearch (isochrone.mass, isochrone.nEntries, zamsMass);

    for (filt = 0; filt < N_CHAB_FILTS; filt++)
    {
        if (useFilt[filt])
        {
            globalMags[filt] = linInterpExtrap (isochrone.mass[m], isochrone.mass[m + 1], isochrone.mag[m][filt], isochrone.mag[m + 1][filt], zamsMass);
            if (fabs (globalMags[filt]) < EPS)
                globalMags[filt] = 999.99;
        }
    }

    return zamsMass;
}

double ChabMsModel::wdPrecLogAge (double FeH, double thisY, double zamsMass)
{

    return 0.0;

    //printf("Chaboyer models unavaible at this time.  Have a nice day!\n (drv_c_AGB_m)");
    exit (1);
}

static void initIso (struct cIsochrone *newIso)
{

    int i, f;

    newIso->age = 0.0;
    newIso->logAge = 0.0;
    newIso->FeH = 0.0;
    newIso->Y = 0.0;
    newIso->AGBt = 0.0;
    newIso->numEeps = 0;
    for (i = 0; i < MAX_CHAB_ENTRIES; i++)
    {
        newIso->eeps[i] = 0;
        newIso->mass[i] = 0.0;
        for (f = 0; f < N_CHAB_FILTS; f++)
            newIso->mag[i][f] = 99.9;
    }
}

static void calcCoeff (double a[], double b[], double x)
{
    b[0] = (a[1] - x) / (a[1] - a[0]);
    b[1] = (x - a[0]) / (a[1] - a[0]);
    return;
}
