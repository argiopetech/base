#include <iostream>
#include <vector>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Cluster.hpp"
#include "Star.hpp"

#include "evolve.hpp"
#include "structures.hpp"

using std::vector;
using std::cerr;
using std::endl;

// Swaps two mass entries in a global isochrone (n and n-1)
void swapGlobalEntries (struct globalIso &thisIso, const vector<int> &filters, int n)
{
    assert(!filters.empty()); // Just in case

    int tempEep;
    double tempMass;
    vector<double> tempMag(filters.back() + 1); // One larger than the largest filter

    tempMass = thisIso.mass.at(n);
    tempEep = thisIso.eep.at(n);
    for (auto f : filters)
        tempMag.at(f) = thisIso.mag.at(n).at(f);

    thisIso.mass.at(n) = thisIso.mass.at(n - 1);
    thisIso.eep.at(n) = thisIso.eep.at(n - 1);
    for (auto f : filters)
            thisIso.mag.at(n).at(f) = thisIso.mag.at(n - 1).at(f);

    thisIso.mass.at(n - 1) = tempMass;
    thisIso.eep.at(n - 1) = tempEep;
    for (auto f : filters)
        thisIso.mag.at(n - 1).at(f) = tempMag.at(f);
}
