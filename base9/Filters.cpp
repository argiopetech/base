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
//      // UBVRIHJK Filters used through December, 2019
//      { "U", 1.569 }, // Cardelli, Clayton, Mathis 1989, table 3
//      { "B", 1.337 }, // yields A_u -> A_k = f(A_v), for standard filters
//      { "V", 1.000 },
//      { "R", 0.751 },
//      { "I", 0.479 },
//      { "J", 0.282 },
//      { "H", 0.190 },
//      { "K", 0.114 },

//      // UBVRIHJK Filters updated January, 2020 from Padova stellar evol webpage
//      { "U", 1.54892 }, // (see http://stev.oapd.inaf.it/cgi-bin/cmd)
//      { "B", 1.29719 },
//      { "V", 1.00600 },
//      { "R", 0.81512 },
//      { "I", 0.60329 },
//      { "J", 0.29553 },
//      { "H", 0.18201 },
//      { "K", 0.11554 },

        // UBVRIHJK Filters updated June, 2022 from Padova stellar evol webpage
        { "U", 1.55814 }, // (see http://stev.oapd.inaf.it/cgi-bin/cmd_3.6)
        { "B", 1.32616 },
        { "V", 1.00096 },
        { "R", 0.80815 },
        { "I", 0.59893 },
        { "J", 0.28688 },
        { "H", 0.18103 },
        { "K", 0.11265 },

        // New DSED 'K' (from prior to Dec, 2019 and still available)
        { "Ks", 0.114 },
        { "Kp", 0.114 },

//      // SDSS Filters used through December, 2019
//      { "u", 5.155 / 3.1 }, // Stoughton et al. (2002, AJ, 123, 485)
//      { "g", 3.793 / 3.1 }, // Table 22, which gives Afilter/E(B-V )
//      { "r", 2.751 / 3.1 }, // We use Afilter/Av, so all are divided
//      { "i", 2.086 / 3.1 }, // by Rv = Av/E(B-V) = 3.1, consistent
//      { "z", 1.479 / 3.1 }, // Cardelli et al. (1989).

//      // SDSS Filters updated January, 2020 from Padova stellar evol webpage
//      { "u", 1.56906 },
//      { "g", 1.20585 },
//      { "r", 0.87122 },
//      { "i", 0.68319 },
//      { "z", 0.49246 },

        // SDSS Filters updated June, 2022 from Padova stellar evol webpage
        { "u", 1.57465 },
        { "g", 1.22651 },
        { "r", 0.86639 },
        { "i", 0.68311 },
        { "z", 0.48245 },

//      // ACS Filters used through December, 2019
//      { "F435W", 4.081 / 3.1 }, // they used R=3.1; also derived via Cardelli et al. values
//      { "F475W", 3.634 / 3.1 }, // Ext. Ratios A(P)/E(B-V) in ACS/WFC Filters for diff. SEDs
//      { "F550M", 3.042 / 3.1 }, // SED F435W F475W F550M F555W F606W F625W F775W F814W
//      { "F555W", 3.177 / 3.1 }, // O5  4.192 3.773 3.052 3.233 2.936 2.673 2.005 1.864
//      { "F606W", 2.809 / 3.1 }, // G2  4.081 3.634 3.042 3.177 2.809 2.637 1.982 1.825
//      { "F625W", 2.637 / 3.1 }, // M0  3.994 3.555 3.030 3.115 2.716 2.616 1.965 1.796
//      { "F775W", 1.982 / 3.1 },
//      { "F814W", 1.825 / 3.1 }, // using values for G2 star, usually ~2% of O5 or M0 value

//      // HST ACS WFC filters updated January, 2020 from Padova stellar evol webpage
//      { "F435W",  1.31785 },
//      { "F475W",  1.19119 },
//      { "F502N",  1.13049 },
//      { "F550M",  0.97994 },
//      { "F555W",  1.03452 },
//      { "F606W",  0.92246 },
//      { "F625W",  0.85085 },
//      { "F658N",  0.80947 },
//      { "F660N",  0.80744 },
//      { "F775W",  0.65375 },
//      { "F814W",  0.60593 },
//      { "F850LP", 0.48440 },

        // HST ACS WFC filters updated June, 2022 from Padova stellar evol webpage
        { "F435W",  1.33879 },
        { "F475W",  1.21179 },
        { "F555W",  1.03065 },
        { "F606W",  0.90328 },
        { "F625W",  0.84550 },
        { "F775W",  0.65239 },
        { "F814W",  0.59696 },

        // UVIS filters (from prior to Dec, 2019 and still available)
        { "UVf275w", 2.094 }, // Based on calculations of Ax/Av by Rachel Wagner-Kaiser
        { "UVf336w", 1.665 }, // from the theoretical values given in the Cardelli paper.
        { "UVf438w", 1.343 },

//      // HST WFC3 filters updated January, 2020 from Padova stellar evol webpage
//      { "WFC218W1", 2.65245 },
//      { "WFC225W1", 2.32794 },
//      { "WFC275W1", 1.94436 },
//      { "WFC336W",  1.65798 },
//      { "WFC390W",  1.42879 },
//      { "WFC438W",  1.33088 },
//      { "WFC475W",  1.18462 },
//      { "WFC555W",  1.04947 },
//      { "WFC606W",  0.92757 },
//      { "WFC625W",  0.86232 },
//      { "WFC775W",  0.66019 },
//      { "WFC814W",  0.61018 },
//      { "WFC105W",  0.37780 },
//      { "WFC110W",  0.33917 },
//      { "WFC125W",  0.29108 },
//      { "WFC140W",  0.24629 },
//      { "WFC160",   0.20433 },

        // HST WFC3 filters updated June, 2022 from Padova stellar evol webpage
        { "WFC218W",  3.08782 },
        { "WFC225W",  2.62940 },
        { "WFC275W",  2.02499 },
        { "WFC336W",  1.67536 },
        { "WFC390W",  1.44380 },
        { "WFC438W",  1.34148 },
        { "WFC475W",  1.20367 },
        { "WFC555W",  1.04089 },
        { "WFC606W",  0.90941 },
        { "WFC625W",  0.85528 },
        { "WFC775W",  0.65905 },
        { "WFC814W",  0.59845 },
        { "WFC105W",  0.36866 },
        { "WFC110W",  0.31708 },
        { "WFC125W",  0.28148 },
        { "WFC140W",  0.23585 },
        { "WFC160W",  0.20177 },

//      // Gaia filters used through December, 2019
//      { "G",    1.976 }, // Note that these appear too high by ~2x
//      { "G_BP", 2.152 },
//      { "G_RP", 1.610 },

//      // GAIA DR2 Maiz filters updated January, 2020 from Padova stellar evol webpage
//      { "G",      0.86105 },
//      { "G_BPbr", 1.06181 },
//      { "G_BPft", 1.07185 },
//      { "G_BP",   1.07185 }, // this filter not in PARSEC, but here for backward compatibility
//      { "G_RP",   0.65069 },

        // Gaia EDR3 filters updated June, 2022 from Padova stellar evol webpage
        { "G",    0.83627 },
        { "G_BP", 1.08337 },
        { "G_RP", 0.63439 },

//      // Pan-STARRS1 Filters (from Elizabeth) used through December, 2019
//      { "g_ps", 1.155 }, // Wang & Chen 2019, astroph 1904.04575v2
//      { "r_ps", 0.843 }, // Table 3
//      { "i_ps", 0.628 },
//      { "z_ps", 0.487 },
//      { "y_ps", 0.395 },

//      // Pan-STARRS1 Filters updated January, 2020 from Padova stellar evol webpage
//      { "g_ps", 1.16529 },
//      { "r_ps", 0.86813 },
//      { "i_ps", 0.67659 },
//      { "z_ps", 0.51743 },
//      { "y_ps", 0.43092 },
//      { "w_ps", 0.85851 },

        // Pan-STARRS1 Filters updated June, 2022 from Padova stellar evol webpage
        { "g_ps", 1.17994 },
        { "r_ps", 0.86190 },
        { "i_ps", 0.67648 },
        { "z_ps", 0.51296 },
        { "y_ps", 0.42905 },
        { "w_ps", 0.83639 }, // this filter is similar to g_ps + r_ps + i_ps

//      // 2MASS filters updated January, 2020 from Padova stellar evol webpage
//      { "J_2M",  0.29434 },
//      { "H_2M",  0.18128 },
//      { "Ks_2M", 0.11838 },

        // 2MASS filters updated June, 2022 from Padova stellar evol webpage
        { "J_2M",  0.28665 },
        { "H_2M",  0.18082 },
        { "Ks_2M", 0.11675 },

//      // UKIDSS filters updated January, 2020 from Padova stellar evol webpage
//      { "UK_Z", 0.49965 },
//      { "UK_Y", 0.38547 },
//      { "UK_J", 0.28887 },
//      { "UK_H", 0.18353 },
//      { "UK_K", 0.11509 },

        // UKIDSS filters updated June, 2022 from Padova stellar evol webpage
        { "UK_Z", 0.49569 },
        { "UK_Y", 0.38408 },
        { "UK_J", 0.28170 },
        { "UK_H", 0.18251 },
        { "UK_K", 0.11281 },

//      // Spitzer filters added January, 2020 from Padova stellar evol webpage
//      { "IRAC3p6", 0.06706 },
//      { "IRAC4p5", 0.05591 },
//      { "IRAC5p8", 0.04948 },
//      { "IRAC8p0", 0.02518 },
//      { "MIPS24",  0.00000 },
//      { "MIPS70",  0.00000 },
//      { "MIPS160", 0.00000 },

        // Spitzer filters added January, 2020 from Padova stellar evol webpage
        { "IRAC3p6", 0.05228 },
        { "IRAC4p5", 0.03574 },
        { "IRAC5p8", 0.02459 },
        { "IRAC8p0", 0.01433 },
        { "MIPS24",  0.00245 },
        { "MIPS70",  0.00041 },
        { "MIPS160", 0.00012 },

//      // WISE filters added January, 2020 from Padova stellar evol webpage
//      { "W1", 0.07134 },
//      { "W2", 0.05511 },
//      { "W3", 0.00220 },
//      { "W4", 0.00000 },

//      // WISE filters added January, 2020 from Wang & Chen (2019, ApJ, 877, 116, Table 3)
//      { "W1", 0.039 }, // (see https://iopscience.iop.org/article/10.3847/1538-4357/ab1c61/pdf)
//      { "W2", 0.026 },
//      { "W3", 0.040 },
//      { "W4", 0.02  }, // No data, so TvH estimate

        // WISE filters added June, 2022 from Padova stellar evol webpage
        { "W1", 0.05688 },
        { "W2", 0.03427 },
        { "W3", 0.00707 },
        { "W4", 0.00274 },

//      // JWST filters added January, 2020 from Padova stellar evol webpage
//      { "JW070W",  0.74602 },
//      { "JW090W",  0.48870 },
//      { "JW115W",  0.32716 },
//      { "JW150W",  0.21465 },
//      { "JW200W",  0.13579 },
//      { "JW277W",  0.08768 },
//      { "JW444W",  0.06708 },
//      { "JW356W",  0.05692 },
//      { "JW150W2", 0.22445 },
//      { "JW322W2", 0.07878 }

        // JWST filters added June, 2022 from Padova stellar evol webpage
        { "JW070W",  0.74534 },
        { "JW090W",  0.47505 },
        { "JW115W",  0.31943 },
        { "JW150W",  0.20942 },
        { "JW200W",  0.13306 },
        { "JW277W",  0.07836 },
        { "JW356W",  0.05191 },
        { "JW444W",  0.03698 },
        { "JW150W2", 0.17131 },
        { "JW322W2", 0.06007 }

    };
