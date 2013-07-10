#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "evolve.hpp"
#include "linInterp.hpp"
#include "binSearch.hpp"
#include "gBaraffeMag.hpp"

using std::array;
using std::ifstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

// Declared in parent program (simCluster or mcmc, or makeCMD)
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



void getBaraffeMags (const vector<int> &filters, array<double, FILTS> &globalMags, double logAge, double mass)
{
    int a, m, i;
    double massMag[2][N_BAR_FILTS];

    const int bDelta = FILTS - N_BAR_FILTS;

    a = binarySearch (barAge, N_BAR_AGES, logAge);
    m = binarySearch (barMass, N_BAR_MASSES, mass);

    //Interpolate in mass first
    for (i = 0; i < 2; i++)
    {
        for (auto f : filters)
        {
            if (f >= bDelta)
            {
                int fCorr = f - bDelta;
                massMag[i][fCorr] = linInterpExtrap (barMass[m], barMass[m + 1], barMag[a + i][m][fCorr], barMag[a + i][m + 1][fCorr], mass);
            }
        }
    }

    //Interpolate in age
    for (auto f : filters)
    {
        if (f >= bDelta)
        {
            int fCorr = f - bDelta;
            globalMags[f] = linInterpExtrap (barAge[a], barAge[a + 1], massMag[0][fCorr], massMag[1][fCorr], logAge);
        }
    }
    for (int filt = 0; filt < (FILTS - N_BAR_FILTS); filt++)
        globalMags[filt] = 99.99;
}
