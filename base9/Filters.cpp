#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "constants.hpp"
#include "Filters.hpp"

using std::cout;
using std::endl;
using std::string;
using std::unordered_map;
using std::vector;

vector<double> Filters::calcAbsCoeffs(const vector<string>& filters)
{
    vector<double> selectAbsCoeffs;

    for (auto f : filters)
    {
        try
        {
            selectAbsCoeffs.push_back(absCoeffs.at(f));
        }
        catch(std::out_of_range &e)
        {
            cout << "Filter \"" << f << "\" does not have an Av coefficient defined" << endl;
            exit(-1);
        }
    }

    return selectAbsCoeffs;
}

const unordered_map<string, double> Filters::absCoeffs =
    {
        // UBVRIHJK Filters
        { "U", 1.569 }, // Cardelli, Clayton, Mathis 1989, table 3
        { "B", 1.337 }, // yields A_u -> A_k = f(A_v), for standard filters
        { "V", 1.000 },
        { "R", 0.751 },
        { "I", 0.479 },
        { "J", 0.282 },
        { "H", 0.190 },
        { "K", 0.114 },

        // New DSED 'K'
        { "Ks", 0.114 },
        { "Kp", 0.114 },

        // SDSS Filters
        { "u", 5.155 / 3.1 }, // Stoughton et al. (2002, AJ, 123, 485)
        { "g", 3.793 / 3.1 }, // Table 22, which gives Afilter/E(B-V )
        { "r", 2.751 / 3.1 }, // We use Afilter/Av, so all are divided
        { "i", 2.086 / 3.1 }, // by Rv = Av/E(B-V) = 3.1, consistent
        { "z", 1.479 / 3.1 }, // Cardelli et al. (1989).

        // ACS Filters
        { "F435W", 4.081 / 3.1 }, // they used R=3.1; also derived via Cardelli et al. values
        { "F475W", 3.634 / 3.1 }, // Ext. Ratios A(P)/E(B-V) in ACS/WFC Filters for diff. SEDs
        { "F550M", 3.042 / 3.1 }, // SED F435W F475W F550M F555W F606W F625W F775W F814W
        { "F555W", 3.177 / 3.1 }, // O5  4.192 3.773 3.052 3.233 2.936 2.673 2.005 1.864
        { "F606W", 2.809 / 3.1 }, // G2  4.081 3.634 3.042 3.177 2.809 2.637 1.982 1.825
        { "F625W", 2.637 / 3.1 }, // M0  3.994 3.555 3.030 3.115 2.716 2.616 1.965 1.796
        { "F775W", 1.982 / 3.1 },
        { "F814W", 1.825 / 3.1 }, // using values for G2 star, usually ~2% of O5 or M0 value

        // UVIS filters
        { "UVf275w", 2.094 }, // Based on calculations of Ax/Av by Rachel Wagner-Kaiser
        { "UVf336w", 1.665 }, // from the theoretical values given in the Cardelli paper.
        { "UVf438w", 1.343 },
    };
