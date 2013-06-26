#include <string>
#include <fstream>
#include <iostream>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "evolve.hpp"
#include "linInterp.hpp"
#include "gBaraffeMag.hpp"

using std::ifstream;
using std::string;
using std::cerr;
using std::endl;

// Declared in parent program (simCluster or mcmc, or makeCMD)
extern int verbose, useFilt[FILTS];
extern double globalMags[FILTS];

static double barAge[N_BAR_AGES];
static double barMass[N_BAR_MASSES];
static double barMag[N_BAR_AGES][N_BAR_MASSES][N_BAR_FILTS];

void loadBaraffe (string path)
{

    int m, a, mStart, filt;
    ifstream pBaraffe;
    char line[1000];
    string tempFile;
    char ind = '0';
    double tempAge;
    float ignore;

    tempFile = path + "COND03_spitzer_0.2-8.2Gyr";

    pBaraffe.open(tempFile);

    if (!pBaraffe)
    {
        cerr << "\n\nFile " << tempFile << " was not found -- Exiting" << endl;
        exit (1);
    }

    //Skip header lines
    a = -1;

    while (pBaraffe.getline(line, 1000))
    {
        sscanf (line, " %c ", &ind);
        if (ind == 't')
        {
            sscanf (line, "%*s %*s %*s %lf ", &tempAge);
            barAge[++a] = log10 (tempAge * 1e9);
            for (m = 0; m < 3; m++)
                pBaraffe.getline(line, 1000);
            if (barAge[a] > log10 (3.6e9))
                mStart = 2;
            else if (barAge[a] > log10 (1.4e9))
                mStart = 1;
            else
                mStart = 0;
            for (m = mStart; m < N_BAR_MASSES; m++)
            {
                pBaraffe >> barMass[m] >> ignore >> ignore >> ignore;
                for (filt = 0; filt < N_BAR_FILTS; filt++)
                    pBaraffe >> barMag[a][m][filt];
            }
        }
    }

    // Kludge to fix missing low mass entries in higher age isochrones
    for (a = 23; a < N_BAR_AGES; a++)
        for (filt = 0; filt < N_BAR_FILTS; filt++)
            barMag[a][1][filt] = barMag[a][2][filt] + (barMag[a - 1][1][filt] - barMag[a - 1][2][filt]);
    for (a = 6; a < N_BAR_AGES; a++)
        for (filt = 0; filt < N_BAR_FILTS; filt++)
            barMag[a][0][filt] = barMag[a][1][filt] + (barMag[a - 1][0][filt] - barMag[a - 1][1][filt]);

    pBaraffe.close();

}



void getBaraffeMags (double logAge, double mass)
{
    int a, m, filt, i;
    double massMag[2][N_BAR_FILTS];

    int binarySearch (double *searchArray, int size, double searchItem);

    a = binarySearch (barAge, N_BAR_AGES, logAge);
    m = binarySearch (barMass, N_BAR_MASSES, mass);

    //Interpolate in mass first
    for (i = 0; i < 2; i++)
    {
        for (filt = 0; filt < N_BAR_FILTS; filt++)
        {
            if (useFilt[filt + 8])
            {
                massMag[i][filt] = linInterpExtrap (barMass[m], barMass[m + 1], barMag[a + i][m][filt], barMag[a + i][m + 1][filt], mass);
            }
        }
    }

    //Interpolate in age
    for (filt = 0; filt < N_BAR_FILTS; filt++)
    {
        if (useFilt[filt + 8])
        {
            globalMags[filt + 8] = linInterpExtrap (barAge[a], barAge[a + 1], massMag[0][filt], massMag[1][filt], logAge);
        }
    }
    for (filt = 0; filt < 8; filt++)
        globalMags[filt] = 99.99;
}
