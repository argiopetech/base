#include <string>
#include <iostream>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "Cluster.hpp"
#include "Star.hpp"

#include "LinearTransform.hpp"
#include "binSearch.hpp"
#include "GirardiMsModel.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;

static int gNumEntries[N_GIR_Z];
static int gBoundary[2][N_GIR_Z][N_GIR_AGES];

static double gFeH    [N_GIR_Z];
static double gIsoMass[N_GIR_Z][MAX_GIR_ENTRIES];
static double gIsoMag [N_GIR_Z][N_GIR_FILTS][MAX_GIR_ENTRIES];
static double gLogAge [N_GIR_Z][N_GIR_AGES];
static double gAGBt   [N_GIR_Z][N_GIR_AGES];

static void calcCoeff (double a[], double b[], double x);

// Set by derive_AGBt_zmass and used by getGirardiMags
static int iAge, iFeH;
static double currentAge, currentFeH;

static char tempFile[100];
static void getFileName (string path, int z, MsFilter filterSet);


void GirardiMsModel::loadModel (string path, MsFilter filterSet)
{

    /****************************************************************************
     ***********************    filterSet = UBVRIJHK      ***********************
     *********************** Girardi models, standard mags **********************
     ** For simplicity, the 8 input files are assumed to have 50 ages each,    **
     ** which is the minimum set among all of these isochrones.  These ages    **
     ** span the range log(age) = 7.80 to 10.25.                               **
     ** The eight files and their metallicities are                            **
     **   iso_stan_z0.50.dat     Z = 0.0                                       **
     **   iso_stan_z0001.50.dat  Z = 0.0001                                    **
     **   iso_stan_z0004.50.dat  Z = 0.0004                                 **
     **   iso_stan_z001.50.dat   Z = 0.001                                  **
     **   iso_stan_z004.50.dat   Z = 0.004                                  **
     **   iso_stan_z008.50.dat   Z = 0.008                                  **
     **   iso_stan_z019.50.dat   Z = 0.019                                  **
     **   iso_stan_z030.50.dat   Z = 0.030                                  **
     ***************************************************************************/

    /****************************************************************************
     **********************        filterSet == ACS       ***********************
     ********************** Girardi models, ACS mags ****************************
     ** For simplicity, the input files are assumed to have 50 ages, just as   **
     ** the above-used Girardi standard mag isochrones.  This meant reducing   **
     ** the 50-74 ages by limiting them to log(age) = 7.80 to 10.25, which is  **
     ** fine for now as the rejected isochrones are younger than any of the    **
     ** clusters for which we have data.  The metallicities are also all the   **
     ** same.  These are the same isochrones, just convolved by Girardi with   **
     ** different filters.  The file names are different only in that the      **
     ** "stan" becomes "acs".                                                  **
     ***************************************************************************/

    int z, a, i;
    double thisLogAge, lastLogAge;
    FILE *pGirardi;
    char line[240];

    if (filterSet != MsFilter::UBVRIJHK && filterSet != MsFilter::SDSS && filterSet != MsFilter::ACS)
    {
        cerr << "\nFilter set " << static_cast<int>(filterSet) << " not available on Girardi models.  Exiting..." << endl;
        exit (1);
    }

    for (z = 0; z < N_GIR_Z; z++)
    {                           // foreach Girardi metallicity/isochrone file

        for (a = 0; a < N_GIR_AGES; a++)
        {                               // initialize age/boundary pointers
            gBoundary[LOW][z][a] = 0;
            gBoundary[HIGH][z][a] = 0;
            gLogAge[z][a] = 0.0;
        }

        for (i = 0; i < MAX_GIR_ENTRIES; i++)
        {                               // initialize array of model parameters
            gIsoMass[z][i] = 0.0; // ZAMS mass
            for (int filt = 0; filt < N_GIR_FILTS; filt++)
                gIsoMag[z][filt][i] = 99.999;   // U through K absolute mags
        }

        getFileName (path, z, filterSet);
        //fscanf(pModelList,"%s",tempFile);                                 // work on one Girardi model at a time
        if ((pGirardi = fopen (tempFile, "r")) == NULL)
        {                               // open file
            cerr << "\nFile " << tempFile << " was not found - exiting" << endl;
            exit (1);
        }


        i = 0;
        a = 0;                  // a = [0,49], gLogAge contains the 50 ages
        lastLogAge = 0.0;
        while (fgets (line, 240, pGirardi) != NULL && i < MAX_GIR_ENTRIES)
        {                               // load first Z=# Girardi model for all ages
            if (line[0] != '#')
            {
                if (filterSet == MsFilter::UBVRIJHK)
                {                       // Girardi UBVRIJHK isocrhones
                    sscanf (line, "%lf %lf %*f %*f %*f %*f %*f %lf %lf %lf %lf %lf %lf %lf %lf %*f", &thisLogAge, &gIsoMass[z][i], &gIsoMag[z][0][i], &gIsoMag[z][1][i], &gIsoMag[z][2][i], &gIsoMag[z][3][i], &gIsoMag[z][4][i], &gIsoMag[z][5][i], &gIsoMag[z][6][i], &gIsoMag[z][7][i]);
                }
                // for filterSet == ACS, use same set of variables but now have F435W F475W F550M F555W F606W F625W F775W F814W
                // absolute mags, instead of UBVRIJHK absolute mags
                if (filterSet == MsFilter::ACS)
                {                       // Girardi hST/ACS/WF isochrones
                    sscanf (line, "%lf %lf %*f %*f %*f %*f %*f %lf %lf %lf %lf %lf %lf %*f %*f %lf %lf %*f %*f %*f", &thisLogAge, &gIsoMass[z][i], &gIsoMag[z][0][i], &gIsoMag[z][1][i], &gIsoMag[z][2][i], &gIsoMag[z][3][i], &gIsoMag[z][4][i], &gIsoMag[z][5][i], &gIsoMag[z][6][i], &gIsoMag[z][7][i]);
                }
                if (fabs (lastLogAge - thisLogAge) > EPS)
                {                       // find model boundaries for ease/speed later
                    gBoundary[LOW][z][a] = i;   // gBoundary contains beginning boundary of age
                    gLogAge[z][a] = thisLogAge;
                    if (a > 0)
                    {
                        gAGBt[z][a - 1] = gIsoMass[z][i - 1];
                        gBoundary[HIGH][z][a - 1] = i - 1;
                    }

                    lastLogAge = thisLogAge;
                    a++;
                }
                i++;
            }
        }
        gAGBt[z][a - 1] = gIsoMass[z][i - 1];     // Add last entry to AGBt table
        gBoundary[HIGH][z][a - 1] = i - 1;
        gNumEntries[z] = i;
    }

    gFeH[0] = -5.0;             // Girardi isochrone metallicities: close enough for Z=0.0
    gFeH[1] = -2.278754;                //                                  equiv to Z=0.0001
    gFeH[2] = -1.676694;                //                                             0.0004
    gFeH[3] = -1.278754;                //                                             0.001
    gFeH[4] = -0.676694;                //                                             0.004
    gFeH[5] = -0.375664;                //                                             0.008
    gFeH[6] = 0.0;              // (i.e., Z_solar=0.019)                       0.019
    gFeH[7] = 0.198368;         //                                             0.030

    ageLimit.first = gLogAge[0][0];
    ageLimit.second = gLogAge[0][N_GIR_AGES - 1];

    fclose (pGirardi);

}

static void getFileName (string path, int z, MsFilter filterSet)
{

    char fileNames[][5] = { "0", "0001", "0004", "001", "004", "008", "019", "030" };

    strcpy (tempFile, "");
    strcat (tempFile, path.c_str());
    if (filterSet == MsFilter::ACS)
        strcat (tempFile, "gIsoACS/iso_acs_z");
    else
        strcat (tempFile, "gIsoStan/iso_stan_z");
    strcat (tempFile, fileNames[z]);
    strcat (tempFile, ".50.dat");

}


double GirardiMsModel::deriveAgbTipMass (const std::vector<int> &filters, double newFeH, double, double newAge)
/****************************************************************************************
last update: 12nov07

Derive AGBt mass (actually the ZAMS mass for the appropriate AGBt star) for a given
cluster age, interpolating in isochrones as necessary.
****************************************************************************************/
{
    double interpAge[2], interpFeH[2];

    if ((newAge < 7.80)
     || (newAge > 10.25)
     || (newFeH < gFeH[0])
     || (newFeH > gFeH[N_GIR_Z - 1]))
    {
        //     log << ("\n Requested FeH too high. (drv_g_AGB_m.c)");
        return 0.0;
    }

    interpAge[HIGH] = ceil (20. * newAge) / 20.;        // round up to nearest 0.05
    interpAge[LOW] = interpAge[HIGH] - 0.05;

    iAge = (int) (rint ((interpAge[LOW] - 7.8) * 20));  // In this line, PFM. Takes the interpAge which we previously rounded to the next lowest 0.05, subtracts the minimum age of the Girardi isochrones, multiplies by the integral of the Girardi logAge step size, and rounds to an integer (which conveniently ends up being the age index).
    iFeH = binarySearch (gFeH, N_GIR_Z, newFeH);        // function returns lower bound

    interpFeH[LOW] = gFeH[iFeH];
    interpFeH[HIGH] = gFeH[iFeH + 1];

    double b[2], d[2];
    calcCoeff (&gFeH[iFeH], d, newFeH);
    calcCoeff (&gLogAge[iFeH][iAge], b, newAge);
/*
    AGBt_zmass_lo = linearTransform<TransformMethod::Interp>(interpAge[LOW], interpAge[HIGH], gAGBt[iFeH][iAge], gAGBt[iFeH][iAge + 1], newAge).val;
    AGBt_zmass_hi = linearTransform<TransformMethod::Interp>(interpAge[LOW], interpAge[HIGH], gAGBt[iFeH + 1][iAge], gAGBt[iFeH + 1][iAge + 1], newAge).val;

    AGBt_zmass = linearTransform<TransformMethod::Interp>(interpFeH[LOW], interpFeH[HIGH], AGBt_zmass_lo, AGBt_zmass_hi, newFeH).val;
*/
    currentAge = newAge;
    currentFeH = newFeH;

    // Discover the range of the lower age isochrone
    int newimax = MAX_GIR_ENTRIES;
    for (int a = iAge; a < iAge + 2; a++)
    {
        for (int z = iFeH; z < iFeH + 2; z++)
        {
            int nMassPoints = gBoundary[HIGH][z][a] - gBoundary[LOW][z][a] + 1;

            if (nMassPoints < newimax)
                newimax = nMassPoints;
        }
    }

    for (int m = 0; m < newimax; m++)
    {
        isochrone.mass[m] = 0.0;
        for (auto f : filters)
            if (f < N_GIR_FILTS)
                isochrone.mag[m][f] = 0.0;

        for (int a = 0; a < 2; a++)
        {
            for (int z = 0; z < 2; z++)
            {
                const int entry = gBoundary[LOW][iFeH + z][iAge + a] + m;
                isochrone.mass[m] += b[a] * d[z] * gIsoMass[iFeH + z][entry];

                for (auto f : filters)
                {
                    if (f < N_GIR_FILTS)
                    {
                        isochrone.mag[m][f] += b[a] * d[z] * gIsoMag[iFeH + z][f][entry];
                    }
                }
            }
        }

        // Sometimes the interpolation process can leave the
        // mass entries out of order.  This swaps them so that
        // the later mass interpolation can function properly
        if (m > 0)
        {
            int n = m;
            while (isochrone.mass[n] < isochrone.mass[n - 1] && n > 0)
            {
                swapGlobalEntries (isochrone, filters, n);
                n--;
            }
        }
    }

    isochrone.nEntries = newimax;
    isochrone.logAge = newAge;
    isochrone.AgbTurnoffMass = isochrone.mass[isochrone.nEntries - 1];

    return isochrone.AgbTurnoffMass;
}


double GirardiMsModel::wdPrecLogAge (double thisFeH, double zamsMass)
/*************************************************************************************
last update: 12nov07

Determine WD precursor age by 2-D interpolating among the AGBt mass versus age values.
Note that the appropriate AGBt mass and lifetime is not the ZAMS mass and lifetime of
the star currently at the AGBt, but rather refers to the properties of the potentially
higher mass and younger AGBt star that was the WD precursor.
*************************************************************************************/
{


    int thisIndexAge[2], f;
    double logAge[2], AGBt_zmass[2], wdPrecLogAge[2], FeHLo, FeHHi, temp;

    AGBt_zmass[HIGH] = AGBt_zmass[LOW] = 0.0;

    FeHLo = gFeH[iFeH];
    FeHHi = gFeH[iFeH + 1];

    // Find wdPrecLogAge for the lower and upper metallicity cases
    for (f = 0; f < 2; f++)
    {

        if (zamsMass < gAGBt[iFeH + f][N_GIR_AGES - 1])
        {                               // possible if cluster older than logAge=9.0
            wdPrecLogAge[f] = gLogAge[iFeH + f][N_GIR_AGES - 1];        // FOR NOW just use logAge = 9.0
            //     log << (" %.3f Mo < smallest AGBt (%.2f) model mass.  Setting precursor log age to %.3f.\n", zamsMass, gAGBt[iFeH + f][N_GIR_AGES - 1], wdPrecLogAge[f]);
        }
        else if (zamsMass > gAGBt[iFeH + f][0])
        {
            wdPrecLogAge[f] = -2.7 * log10 (zamsMass / gAGBt[iFeH + f][0]) + gLogAge[iFeH + f][0];

            //     log << (" %.3f Mo > largest AGBt (%.2f) model mass.  Extrapolating precursor log age.\n", zamsMass, gAGBt[iFeH + f][0]);
        }

        else
        {
            thisIndexAge[LOW] = reverseBinarySearch (gAGBt[iFeH + f], N_GIR_AGES, zamsMass);    // Because masses are in reverse order
            thisIndexAge[HIGH] = thisIndexAge[LOW] + 1;

            logAge[LOW] = gLogAge[iFeH + f][thisIndexAge[LOW]];
            logAge[HIGH] = gLogAge[iFeH + f][thisIndexAge[HIGH]];

            AGBt_zmass[LOW] = gAGBt[iFeH + f][thisIndexAge[LOW]];
            AGBt_zmass[HIGH] = gAGBt[iFeH + f][thisIndexAge[HIGH]];

            // Linearly interpolate in mass
            wdPrecLogAge[f] = linearTransform<TransformMethod::Interp>(AGBt_zmass[LOW], AGBt_zmass[HIGH], logAge[LOW], logAge[HIGH], zamsMass).val;
        }
    }

    // Linearly interpolate in FeH
    temp = linearTransform<TransformMethod::Interp>(FeHLo, FeHHi, wdPrecLogAge[LOW], wdPrecLogAge[HIGH], thisFeH).val;

    return temp;

}

// a and b are 1-d arrays with two elements
// a contains the two bounding values to be interpolated
// x is the value to be interpolated at
// b returns the coefficients
static void calcCoeff (double a[], double b[], double x)
{
    b[0] = (a[1] - x) / (a[1] - a[0]);
    b[1] = (x - a[0]) / (a[1] - a[0]);
    return;
}
