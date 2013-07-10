#include <array>
#include <iostream>
#include <string>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <gsl/gsl_errno.h>      // need these two lines for gnu interpolation -- TvH
#include <gsl/gsl_spline.h>

#include "evolve.hpp"
#include "DsedMsModel.hpp"
#include "binSearch.hpp"
#include "linInterp.hpp"
#include "Matrix.hpp"

using std::array;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

extern double globalMags[FILTS];
extern struct globalIso isochrone;

// defined in evolve.c, used to determine the normalization
// for the prior on age.
extern double ageLimit[2];

struct dIsochrone
{
    dIsochrone()
    {
        eeps.fill(0);
        mass.fill(0);

        for (auto a : mag)
            a.fill(99.9);

//        assert(false);
    }

    Matrix<double, MAX_DSED_ENTRIES, N_DSED_FILTS> mag;

    array<int, MAX_DSED_ENTRIES> eeps;

    array<double, MAX_DSED_ENTRIES> mass;

    int numEeps = 0;

    double age = 0.0; //In Gyr
    double logAge = 0.0;
    double FeH = 0.0;
    double AGBt;
};

static int iFeH, iAge;
static double dFeH[N_DSED_Z];
static double dLogAge[N_DSED_Z][N_DSED_AGES], dAge[N_DSED_Z][N_DSED_AGES];
static struct dIsochrone dIso[N_DSED_Z][N_DSED_AGES];
static double dAGBt[N_DSED_Z][N_DSED_AGES];
static double dAgbCoeff[N_DSED_Z][2] = { {-3.3734996, 0.8357878},
                                         {-3.3849670, 0.9768837},
                                         {-3.393597, 1.076594},
                                         {-3.506194, 1.280817},
                                         {-3.520550, 1.139945},
                                         {-3.647955, 1.039043},
                                         {-3.6150796, 0.8990882},
                                         {-3.6248249, 0.8607815},
                                         {-3.5500060, 0.6773844}
};

static void calcCoeff (double a[], double b[], double x);

void DsedMsModel::loadModel (string path, MsFilter filterSet)
{

    FILE *pDsed;                        // = NULL;
    int z, a, i, f;
    char line[240];
    string tempFile;

    if (filterSet != MsFilter::SDSS && filterSet != MsFilter::UBVRIJHK)
    {
        cerr << "\nFilter set " << static_cast<int>(filterSet) << " not available on DSED models.  Exiting..." << endl;
        exit (1);
    }


    for (z = 0; z < N_DSED_Z; z++)
    {                           // foreach Dsed metallicity/isochrone file
        dFeH[z] = 0.0;

        for (a = 0; a < N_DSED_AGES; a++)
        {                               // initialize age/boundary pointers
            dLogAge[z][a] = 0.0;
            dAge[z][a] = 0.0;
        }
        a = -1;                 // a = [0,18], dLogAge contains the 52 ages
        i = 0;

        for (f = 0; f < 2; f++)
        {                               // Each metallicity has 2 files, 0.25-1 Gyr and 1-15 Gyr
            tempFile = getFileName (path, z, f, filterSet);        // work on one Dsed model at a time
            if ((pDsed = fopen (tempFile.c_str(), "r")) == NULL)
            {                           // open file
                cerr << "\n\n file " << tempFile << " was not found - exiting" << endl;
                exit (1);
            }
            while (fgets (line, 240, pDsed) != NULL)
            {                           // load each Z=# Dsed model for all ages
                if (line[1] == 'M')
                {
                    if (f == 0)
                    {
                        fgets (line, 240, pDsed);
                        sscanf (line, "%*s %*f %*f %*f %*f %lf", &dFeH[z]);
                    }
                }
                else if (line[1] == 'A')
                {
                    a++;
                    // Skip the first entry of the second file, since it
                    // duplicates the last entry of the first file
                    if (a == 16)
                    {
                        while (line[0] != '\n')
                            fgets (line, 240, pDsed);
                        fgets (line, 240, pDsed);
                        fgets (line, 240, pDsed);
                    }
                    sscanf (line, "#AGE= %lf EEPS= %d", &(dIso[z][a].age), &(dIso[z][a].numEeps));
                    dLogAge[z][a] = dIso[z][a].logAge = log10 (dIso[z][a].age * (1e9));
                    dAge[z][a] = dIso[z][a].age;
                    dIso[z][a].FeH = dFeH[z];
                    i = 0;
                }
                else if (line[0] != '#' && line[0] != '\n')
                {
                    sscanf (line, "%d %lf %*f %*f %*f %lf %lf %lf %lf %lf %lf %lf %lf", &dIso[z][a].eeps.at(i), &dIso[z][a].mass.at(i), &dIso[z][a].mag[i][0], &dIso[z][a].mag[i][1], &dIso[z][a].mag[i][2], &dIso[z][a].mag[i][3], &dIso[z][a].mag[i][4], &dIso[z][a].mag[i][5], &dIso[z][a].mag[i][6], &dIso[z][a].mag[i][7]);

                    if (i == dIso[z][a].numEeps - 1)
                    {
                        dIso[z][a].AGBt = dIso[z][a].mass.at(i);
                    }
                    i++;
                }
            }
        }

        // Lop off the extra bit of evolution in the lower age models
        for (a = 0; a < N_DSED_AGES; a++)
        {
            if (dIso[z][a].numEeps > 280)
                dIso[z][a].numEeps -= dIso[z][a].eeps.at(dIso[z][a].numEeps - 1) - 220;
            dAGBt[z][a] = dIso[z][a].mass.at(dIso[z][a].numEeps - 1);
        }
    }

    ageLimit[0] = dLogAge[0][0];
    ageLimit[1] = dLogAge[0][N_DSED_AGES - 1];

    //Load in JHK from the UBVRIJHK models
    if (filterSet == MsFilter::SDSS)
    {
        filterSet = MsFilter::UBVRIJHK;

        for (z = 0; z < N_DSED_Z; z++)
        {                               // foreach Dsed metallicity/isochrone file
            a = -1;                     // a = [0,18], dLogAge contains the 52 ages
            i = 0;

            for (f = 0; f < 2; f++)
            {                           // Each metallicity has 2 files, 0.25-1 Gyr and 1-15 Gyr
                getFileName (path, z, f, filterSet);    // work on one Dsed model at a time

                if ((pDsed = fopen (tempFile.c_str(), "r")) == NULL)
                {                       // open file
                    cerr << "\n\n file " << tempFile << " was not found - exiting" << endl;
                    exit (1);
                }

                while (fgets (line, 240, pDsed) != NULL)
                {                       // load each Z=# Dsed model for all ages
                    if (line[1] == 'A')
                    {
                        a++;
                        // Skip the first entry of the second file, since it
                        // duplicates the last entry of the first file

                        if (a == 16)
                        {
                            while (line[0] != '\n')
                                fgets (line, 240, pDsed);
                            fgets (line, 240, pDsed);
                            fgets (line, 240, pDsed);
                        }
                        i = 0;
                    }
                    else if (line[0] != '#' && line[0] != '\n')
                    {
                        sscanf (line, "%*d %*f %*f %*f %*f %*f %*f %*f %*f %*f %lf %lf %lf", &dIso[z][a].mag[i][5], &dIso[z][a].mag[i][6], &dIso[z][a].mag[i][7]);
                        i++;
                    }
                }
            }
        }
    }
}


string DsedMsModel::getFileName (string path, int z, int f, MsFilter filterSet)
{

    const array<string, 9> fileNames = {{ "m25", "m20", "m15", "m10", "m05", "p00", "p02", "p03", "p05" }};

    if (filterSet == MsFilter::SDSS)
        path += "sdss/feh";
    else
        path += "jc2mass/feh";

    path += fileNames[z];
    path += "afep0.";

    if (filterSet == MsFilter::SDSS)
        path += "ugriz";
    else
        path += "jc2mass";

    if (!f)
        path += "_2";

    return path;
}

// Interpolates between isochrones for two ages using linear interpolation
// Must run loadDsed() first for this to work.
// Currently ignores newY
double DsedMsModel::deriveAgbTipMass (const std::vector<int> &filters, double newFeH, double newY, double newLogAge)
{

    int newimax = 500, newimin = 0, ioff[2][2], neweep;
    int z = 0, a = 0, m = 0, n = 0;
    double newAge = exp10 (newLogAge) / 1e9;
    double b[2], d[2];

    iAge = -1;
    iFeH = -1;

    if (newLogAge < dLogAge[0][0])
    {
        //     log << ("\n Requested age (%.3f) too young. (gDsedMag.c)", newLogAge);
        return 0.0;
    }

    if (newLogAge > dLogAge[N_DSED_Z - 1][N_DSED_AGES - 1])
    {
        //     log << ("\n Requested age (%.3f) too old. (gDsedMag.c)", newLogAge);
        return 0.0;
    }

    if (newFeH < dFeH[0])
    {
        //     log << ("\n Requested FeH (%.3f) too low. (gDsedMag.c)", newFeH);
        return 0.0;
    }

    if (newFeH > dFeH[N_DSED_Z - 1])
    {
        //     log << ("\n Requested FeH (%.3f) too high. (gDsedMag.c)", newFeH);
        return 0.0;
    }

    // Find the values for each parameter that we will be interpolating
    // between and calculate the interpolation coefficients.
    iFeH = binarySearch (dFeH, N_DSED_Z, newFeH);
    iAge = binarySearch (dAge[iFeH], N_DSED_AGES, newAge);

    calcCoeff (&dFeH[iFeH], d, newFeH);
    calcCoeff (&dAge[iFeH][iAge], b, newAge);

    // Find the minimum and maximum eep values and set a global
    // min and max that is the union of the min and max for each isochrone
    for (a = iAge; a < iAge + 2; a++)
    {
        for (z = iFeH; z < iFeH + 2; z++)
        {
            if (dIso[z][a].eeps.at(0) > newimin)
                newimin = dIso[z][a].eeps.at(0);
            if (dIso[z][a].eeps.at(dIso[z][a].numEeps - 1) < newimax)
                newimax = dIso[z][a].eeps.at(dIso[z][a].numEeps - 1);
        }
    }

    neweep = newimax - newimin + 1;     // = the # of eep points that will be in the new isochrone

    // For each isochrone, find the amount by which
    // the min eep # for that isochrone is
    // offset from the global minimum eep #
    for (a = 0; a < 2; a++)
    {
        for (z = 0; z < 2; z++)
        {
            ioff[z][a] = newimin - dIso[z + iFeH][a + iAge].eeps.at(0);
        }
    }

    // Now for each entry in the new isochrone
    // use the coefficients to calculate the mass
    // and each color at that eep
    for (m = 0; m < neweep; m++)
    {
        isochrone.mass[m] = 0.0;
        for (auto f : filters)
            if (f < N_DSED_FILTS)
                isochrone.mag[m][f] = 0.0;

        for (a = 0; a < 2; a++)
        {
            for (z = 0; z < 2; z++)
            {
                isochrone.mass[m] += b[a] * d[z] * dIso[iFeH + z][iAge + a].mass.at(m + ioff[z][a]);
                for (auto f : filters)
                {
                    if (f < N_DSED_FILTS)
                    {
                        isochrone.mag[m][f] += b[a] * d[z] * dIso[iFeH + z][iAge + a].mag[m + ioff[z][a]][f];
                    }
                }
            }
        }

        isochrone.eep[m] = dIso[iFeH][iAge].eeps.at(m + ioff[0][0]);

        // Sometimes the interpolation process can leave the
        // mass entries out of order.  This swaps them so that
        // the later mass interpolation can function properly
        if (m > 0)
        {
            n = m;
            while (isochrone.mass[n] < isochrone.mass[n - 1] && n > 0)
            {
                swapGlobalEntries (isochrone, filters, n);
                n--;
            }
        }
    }

    //Transfer to globalIsochrone structure
    isochrone.nEntries = neweep;
    isochrone.logAge = newLogAge;
    isochrone.FeH = newFeH;
    isochrone.AgbTurnoffMass = isochrone.mass[isochrone.nEntries - 1];

    return isochrone.AgbTurnoffMass;
}


// Calculates magnitudes for a given mass.
// Must run loadDsed() and deriveDsedAgbTip()
// to load and interpolate an isochrone before this subroutine will work
// Stores output values in external variable globalMags[]
double DsedMsModel::msRgbEvol (const vector<int> &filters, double zamsMass)
{
    int m;

    m = binarySearch (isochrone.mass.data(), isochrone.nEntries, zamsMass);

    for (auto f : filters)
    {
        if (f < N_DSED_FILTS)
        {
            globalMags[f] = linInterpExtrap (isochrone.mass[m], isochrone.mass[m + 1], isochrone.mag[m][f], isochrone.mag[m + 1][f], zamsMass);
            if (fabs (globalMags[f]) < EPS)
                globalMags[f] = 999.99;
        }
    }

    return zamsMass;
}


// Calculates the precursor age for a given wd precursor mass
// Must run loadDsed() and deriveDsedAgbTip()
// to load and interpolate an isochrone before this subroutine will work
double DsedMsModel::wdPrecLogAge (double thisFeH, double thisY, double zamsMass)
{
    int thisIndexAge[2], f;
    double logAge[2], AgbTurnoffMass[2], wdPrecLogAge[2], FeH[2], temp;

    int reverseBinarySearch (double *searchArray, int size, double searchItem);

    AgbTurnoffMass[1] = AgbTurnoffMass[0] = 0.0;

    FeH[0] = dFeH[iFeH];
    FeH[1] = dFeH[iFeH + 1];

    // Find wdPrecLogAge for the lower and upper metallicity cases
    for (f = 0; f < 2; f++)
    {
        if (zamsMass < dAGBt[iFeH + f][N_DSED_AGES - 1]) // possible if cluster older than logAge=9.0
        {
            wdPrecLogAge[f] = dLogAge[iFeH + f][N_DSED_AGES - 1];       // FOR NOW just use logAge = 9.0
            //     log << (" %.3f Mo < smallest AGBt (%.2f) model mass.  Setting precursor log age to %.3f.\n", zamsMass, dAGBt[iFeH + f][N_DSED_AGES - 1], wdPrecLogAge[f]);
        }
        else if (zamsMass > dAGBt[iFeH + f][0])
        {
            wdPrecLogAge[f] = dAgbCoeff[iFeH + f][1] * (pow (log10 (zamsMass), 2) - pow (log10 (dAGBt[iFeH + f][0]), 2)) + dAgbCoeff[iFeH + f][0] * (log10 (zamsMass / dAGBt[iFeH + f][0])) + dLogAge[iFeH + f][0];

            //     log << (" %.3f Mo > largest AGBt (%.2f) model mass.  Extrapolating precursor log age.\n", zamsMass, dAGBt[iFeH + f][0]);
        }
        else
        {

            thisIndexAge[0] = reverseBinarySearch (dAGBt[iFeH + f], N_DSED_AGES, zamsMass);     // Because masses are in reverse order
            thisIndexAge[1] = thisIndexAge[0] + 1;

            logAge[0] = dLogAge[iFeH + f][thisIndexAge[0]];
            logAge[1] = dLogAge[iFeH + f][thisIndexAge[1]];

            AgbTurnoffMass[0] = dAGBt[iFeH + f][thisIndexAge[0]];
            AgbTurnoffMass[1] = dAGBt[iFeH + f][thisIndexAge[1]];

            // Linearly interpolate in mass
            wdPrecLogAge[f] = linInterp (AgbTurnoffMass[0], AgbTurnoffMass[1], logAge[0], logAge[1], zamsMass);
        }
    }

    // Linearly interpolate in FeH
    temp = linInterp (FeH[0], FeH[1], wdPrecLogAge[0], wdPrecLogAge[1], thisFeH);

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
