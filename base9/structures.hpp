#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP

#include <vector>

#include <cassert>

#include "constants.hpp"

struct globalIso
{
    globalIso()
        : globalIso(370, FILTS)
    {;}

    globalIso(int nEntries, int nFilts)
        : nEntries(nEntries)
    {
        assert(nEntries != 0 && nFilts != 0);
    
        mass.resize(nEntries, 0.0);
        eep.resize(nEntries, 0);

        mag.assign(nEntries, std::vector<double>(nFilts, 0.0));
    }

    std::vector< std::vector<double> > mag;

    std::vector<int> eep;
    std::vector<double> mass;

    unsigned int nEntries;

    double AgbTurnoffMass = 0.0;
    double FeH = 0.0;
    double logAge = 0.0;
    double age = 0.0;
    double Y = 0.0;
    double z = 0.0;
};

void swapGlobalEntries (struct globalIso &thisIso, const std::vector<int>&, int n);

#endif
