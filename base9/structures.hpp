#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP

#include <vector>

#include <cassert>

#include "constants.hpp"
#include "Chain.hpp"

struct globalIso
{
    globalIso()
        : globalIso(370, FILTS)
    {;}

    globalIso(int nEntries, int nFilts)
        : nEntries(nEntries)
    {
        assert(nEntries != 0 && nFilts != 0);
    
        mass.resize(nEntries);
        eep.resize(nEntries);

        mag.assign(nEntries, std::vector<double>(nFilts));
    }

    std::vector< std::vector<double> > mag;

    std::vector<int> eep;
    std::vector<double> mass;

    int nEntries;

    double AgbTurnoffMass;
    double FeH;
    double logAge;
    double age;
    double Y;
    double z;
};

enum class blockType
{
    INITIAL_WANDER_BLOCK,
    MASS_RATIO_BLOCK,
    AGE_BLOCK,
    MODULUS_BLOCK,
    AGE_MODULUS_BLOCK,
    FEH_BLOCK,
    Y_BLOCK,
    TURN_FS_SAMPLING_OFF_BLOCK,
    FIELD_STAR_PRIOR_BLOCK,
    TURN_FS_SAMPLING_ON_BLOCK,
    FINAL_WANDER_BLOCK,
    END_BURN_IN_BLOCK,
    POST_BURN_IN_BLOCK
};

void swapGlobalEntries (struct globalIso &thisIso, const std::vector<int>&, int n);

#endif
