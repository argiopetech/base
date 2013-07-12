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
    
        mass.resize(nEntries, 0.0);
        eep.resize(nEntries, 0);

        mag.assign(nEntries, std::vector<double>(nFilts, 0.0));
    }

    std::vector< std::vector<double> > mag;

    std::vector<int> eep;
    std::vector<double> mass;

    uint nEntries;

    double AgbTurnoffMass = 0.0;
    double FeH = 0.0;
    double logAge = 0.0;
    double age = 0.0;
    double Y = 0.0;
    double z = 0.0;
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
