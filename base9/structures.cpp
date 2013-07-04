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
