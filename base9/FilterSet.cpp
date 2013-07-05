#include <array>
#include <iostream>

#include <cstring>

#include "constants.hpp"

using std::array;
using std::cerr;
using std::endl;

array<double, FILTS> calcAbsCoeffs (MsFilterSet filterSet)
{
    if (filterSet == MsFilterSet::UBVRIJHK)
    {
        return { 1.569  // Cardelli, Clayton, Mathis 1989, table 3
                     , 1.337  // yields A_u -> A_k = f(A_v), for standard filters
                     , 1.0
                     , 0.751
                     , 0.479
                     , 0.282
                     , 0.190
                     , 0.114};
    }
    else if (filterSet == MsFilterSet::SDSS)
    {
        return { 5.155 / 3.1    // Stoughton et al. (2002, AJ, 123, 485)
                     , 3.793 / 3.1    // Table 22, which gives Afilter/E(B-V )
                     , 2.751 / 3.1    // We use Afilter/Av, so all are divided
                     , 2.086 / 3.1    // by Rv = Av/E(B-V) = 3.1, consistent
                     , 1.479 / 3.1    // Cardelli et al. (1989).
                     , 0.282  // JHK come from Cardelli (see above)
                     , 0.190
                     , 0.114};
    }
    else if (filterSet == MsFilterSet::ACS)
    {                           // from Table 14, Sirianni et al. (2005, PASP, 117, 1049)
        return { 4.081 / 3.1    // they used R=3.1; also derived via Cardelli et al. values
                     , 3.634 / 3.1    // Ext. Ratios A(P)/E(B-V) in ACS/WFC Filters for diff. SEDs
                     , 3.042 / 3.1    // SED F435W F475W F550M F555W F606W F625W F775W F814W
                     , 3.177 / 3.1    // O5  4.192 3.773 3.052 3.233 2.936 2.673 2.005 1.864
                     , 2.809 / 3.1    // G2  4.081 3.634 3.042 3.177 2.809 2.637 1.982 1.825
                     , 2.637 / 3.1    // M0  3.994 3.555 3.030 3.115 2.716 2.616 1.965 1.796
                     , 1.982 / 3.1
                     , 1.825 / 3.1};    // using values for G2 star, usually ~2% of O5 or M0 value
    }
    else
    {
        cerr << "filterSet " << static_cast<int>(filterSet) << "not found.  Exiting." << endl;
        exit (1);
    }
}


char filterNames[FILTS][10];

void setFilterNames (MsFilterSet filterSet)
{
    int f;

    for (f = 0; f < FILTS; f++)
        strcpy (filterNames[f], "\0");
    if (filterSet == MsFilterSet::UBVRIJHK)
    {
        strcat (filterNames[0], "U");
        strcat (filterNames[1], "B");
        strcat (filterNames[2], "V");
        strcat (filterNames[3], "R");
        strcat (filterNames[4], "I");
        strcat (filterNames[5], "J");
        strcat (filterNames[6], "H");
        strcat (filterNames[7], "K");
    }
    else if (filterSet == MsFilterSet::ACS)
    {
        strcat (filterNames[0], "F435W");
        strcat (filterNames[1], "F475W");
        strcat (filterNames[2], "F550M");
        strcat (filterNames[3], "F555W");
        strcat (filterNames[4], "F606W");
        strcat (filterNames[5], "F625W");
        strcat (filterNames[6], "F775W");
        strcat (filterNames[7], "F814W");
    }
    else if (filterSet == MsFilterSet::SDSS)
    {
        strcat (filterNames[0], "u");
        strcat (filterNames[1], "g");
        strcat (filterNames[2], "r");
        strcat (filterNames[3], "i");
        strcat (filterNames[4], "z");
        strcat (filterNames[5], "J");
        strcat (filterNames[6], "H");
        strcat (filterNames[7], "K");
    }
    else
    {
        cerr << "\nFilterSet " << static_cast<int>(filterSet) << " not available.  Exiting." << endl;
        exit (1);
    }

    strcat (filterNames[8], "IrB");
    strcat (filterNames[9], "IrR");
    strcat (filterNames[10], "SpB1");
    strcat (filterNames[11], "SpB2");
    strcat (filterNames[12], "SpB3");
    strcat (filterNames[13], "SpB4");
}

char *getFilterName (int index)
{
    return filterNames[index];
}
