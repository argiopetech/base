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
void swapGlobalEntries (struct globalIso &thisIso, const vector<int> &filters, int n)
{
    int tempEep;
    double tempMass;
    vector<double> tempMag(thisIso.nFilts);

    tempMass = thisIso.mass[n];
    tempEep = thisIso.eep[n];
    for (auto f : filters)
        tempMag[f] = thisIso.mag[n][f];

    thisIso.mass[n] = thisIso.mass[n - 1];
    thisIso.eep[n] = thisIso.eep[n - 1];
    for (auto f : filters)
            thisIso.mag[n][f] = thisIso.mag[n - 1][f];

    thisIso.mass[n - 1] = tempMass;
    thisIso.eep[n - 1] = tempEep;
    for (auto f : filters)
        thisIso.mag[n - 1][f] = tempMag[f];
}
