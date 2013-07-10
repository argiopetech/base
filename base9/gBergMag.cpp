/********************************************************************************
last update: 19jul10

For a given H or He atmosphere WD Teff and logg, determine U through K (or u
through z) mags from 2-D interpolating in the Bergeron et al. (1995, PASP, 107,
1047) H and He model atmosphere tables.  Interpolating in log(Teff) since
luminosities and colors are all in log form and since the log-log plots are
significantly more linear.  (Note that for all but the hottest WDs the Bergeron
et al. atmospheres span a wider mass range than the Wood model WDs.)

The appropriate magnitudes are put in globalMags[][].
********************************************************************************/

#include <array>
#include <string>
#include <iostream>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "evolve.hpp"
#include "linInterp.hpp"
#include "gBergMag.hpp"
#include "binSearch.hpp"

using std::array;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

// Arrays to hold the models
static double bLogG[BERG_N_DA_LOG_G];
static double bMag[BERG_N_DA_TEFF][BERG_N_DA_LOG_G][BERG_NFILTS];
static double bLogTeff[BERG_N_DA_TEFF];

static double bHeLogG[BERG_N_DB_LOG_G];
static double bHeMag[BERG_N_DB_TEFF][BERG_N_DB_LOG_G][BERG_NFILTS];
static double bHeLogTeff[BERG_N_DB_TEFF];

void loadBergeron (string path, MsFilter filterSet)
{
    int l, t, f;
    FILE *pBergeron;
    char line[1000], tempFile[100] = "\0";
    double tempTeff;

    if (filterSet != MsFilter::SDSS && filterSet != MsFilter::UBVRIJHK)
    {
        cerr << "\nFilter set " << static_cast<int>(filterSet) << " not available on Bergeron models.  Exiting..." << endl;
        exit (1);
    }

    // Open the appropriate file for each mass
    //fscanf(pModelList,"%s",tempFile);
    strcat (tempFile, path.c_str());
    strcat (tempFile, "bergeron/Table_DA.txt");

    if ((pBergeron = fopen (tempFile, "r")) == NULL)
    {
        cerr << "\nFile " << tempFile << " was not found - exiting" << endl;
        exit (1);
    }

    //Skip header lines
    fgets (line, 1000, pBergeron);
    fgets (line, 1000, pBergeron);

    for (l = 0; l < BERG_N_DA_LOG_G; l++)
    {
        for (t = 0; t < BERG_N_DA_TEFF; t++)
        {
            // Teff logg M/Mo Mbol BC U B V R I J H K u g r i z y b-y u-b v-y V-I G-R U-V U-G B-V Age
            fscanf (pBergeron, "%lf %lf %*f %*f %*f ", &tempTeff, &bLogG[l]);
            for (f = 0; f < BERG_NFILTS; f++)
                fscanf (pBergeron, "%lf ", &bMag[t][l][f]);
            if (filterSet == MsFilter::SDSS)
                for (f = 0; f < 5; f++)
                    fscanf (pBergeron, "%lf ", &bMag[t][l][f]);
            fgets (line, 1000, pBergeron);
            bLogTeff[t] = log10 (tempTeff);
        }
    }

    //Still need to add He models
    // Open the appropriate file for each mass
    strcpy (tempFile, "\0");
    strcat (tempFile, path.c_str());
    strcat (tempFile, "bergeron/Table_DB.txt");

    if ((pBergeron = fopen (tempFile, "r")) == NULL)
    {
        cerr << "\n file " << tempFile << " was not found - exiting" << endl;
        exit (1);
    }

    //Skip header lines
    fgets (line, 1000, pBergeron);
    fgets (line, 1000, pBergeron);

    for (l = 0; l < BERG_N_DB_LOG_G; l++)
    {
        for (t = 0; t < BERG_N_DB_TEFF; t++)
        {
            // Teff  log g  M/Mo  Mbol     BC    U      B      V      R      I      J      H      K      u      g      r      i      z      y      b-y    u-b    v-y    V-I    G-R    U-V    U-G    B-V     Age
            fscanf (pBergeron, "%lf %lf %*f %*f %*f ", &tempTeff, &bHeLogG[l]);
            for (f = 0; f < BERG_NFILTS; f++)
                fscanf (pBergeron, "%lf ", &bHeMag[t][l][f]);
            if (filterSet == MsFilter::SDSS)
                for (f = 0; f < 5; f++)
                    fscanf (pBergeron, "%lf ", &bHeMag[t][l][f]);
            bHeLogTeff[t] = log10 (tempTeff);
            fgets (line, 1000, pBergeron);
        }
    }
    fclose (pBergeron);
}


void bergeronTeffToMags (const vector<int> &filters, array<double, FILTS> &globalMags, double wdLogTeff, double wdLogG, int wdType)
{
    int l, t, i;
    double logGMag[2][BERG_NFILTS];

    if (wdType == DA)
    {

        if (wdLogG < bLogG[0])
            l = 0;
        else if (wdLogG >= bLogG[BERG_N_DA_LOG_G - 2])
            l = BERG_N_DA_LOG_G - 2;
        else
            l = (int) floor ((wdLogG - bLogG[0]) / (bLogG[1] - bLogG[0]));
        t = binarySearch (bLogTeff, BERG_N_DA_TEFF, wdLogTeff);

        //Interpolate in logTeff
        for (i = 0; i < 2; i++)
        {
            for (auto f : filters)
            {
                if (f < BERG_NFILTS)
                {
                    logGMag[i][f] = linInterpExtrap (bLogTeff[t], bLogTeff[t + 1], bMag[t][l + i][f], bMag[t + 1][l + i][f], wdLogTeff);
                }
            }
        }

        //Interpolate in log(g)
        for (auto f : filters)
        {
            if (f < BERG_NFILTS)
            {
                globalMags[f] = linInterpExtrap (bLogG[l], bLogG[l + 1], logGMag[0][f], logGMag[1][f], wdLogG);
            }
        }
    }

    else if (wdType == DB)
    {
        if (wdLogG < bHeLogG[0])
            l = 0;
        else if (wdLogG >= bHeLogG[BERG_N_DB_LOG_G - 2])
            l = BERG_N_DB_LOG_G - 2;
        else
            l = (int) floor ((wdLogG - bHeLogG[0]) / (bHeLogG[1] - bHeLogG[0]));
        t = binarySearch (bLogTeff, BERG_N_DB_TEFF, wdLogTeff);

        //Interpolate in logTeff
        for (i = 0; i < 2; i++)
        {
           for (auto f : filters)
            {
                if (f < BERG_NFILTS)
                {
                    logGMag[i][f] = linInterpExtrap (bHeLogTeff[t], bHeLogTeff[t + 1], bHeMag[t][l + i][f], bHeMag[t + 1][l + i][f], wdLogTeff);
                }
            }
        }

        //Interpolate in log(g)
        for (auto f : filters)
        {
            if (f < BERG_NFILTS)
            {
                globalMags[f] = linInterpExtrap (bHeLogG[l], bHeLogG[l + 1], logGMag[0][f], logGMag[1][f], wdLogG);
            }
        }
    }
    else
    {
        // NS EDIT: it's not very nice to exit...
        cerr << "\nInvalid WD type.  Must be 0 (DA) or 1 (DB). Not exiting..." << endl;
        // exit(1);
        for (auto f : filters)
        {
            if (f < BERG_NFILTS)
            {
                globalMags[f] = 99.99;
            }
        }
    }
}
