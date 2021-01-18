#ifndef BERG2020ATMOS_HPP
#define BERG2020ATMOS_HPP

#include <array>
#include <map>
#include <string>
#include <utility>

#include "constants.hpp"
#include "WdAtmosphereModel.hpp"

class Bergeron2020AtmosphereModel : public BergeronAtmosphereModel
{
  public:
    Bergeron2020AtmosphereModel()
    {
        dirName = "bergeron_2020/";

        availableFilters = { "U", "B", "V", "R", "I",                // Johnson-Kron-Cousins UBVRI
                             "J", "H", "Ks",                         // 2MASS JHKs, where Ks means "K short"
                             "Y_MKO", "J_MKO", "H_MKO", "K_MKO",     // Mauna Kea Observatory (MKO) YJHK
                             "W1", "W2", "W3", "W4",                 // Wide-field Infrared Survey Explorer (WISE) W1 to W4
                             "S3.6", "S4.5", "S5.8", "S8.0",         // Spitzer IRAC 3.6, 4.5, 5.8, 8.0 microns
                             "u", "g", "r", "i", "z",                // SDSS
                             "g_ps", "r_ps", "i_ps", "z_ps", "y_ps", // Pan-STARRS
                             "G2", "G2_BP", "G2_RP",                 // Gaia 2
                             "G3", "G3_BP", "G3_RP",                 // Gaia 3
                             "FUV", "NUV" };                         // Galex

        files = {
            {"Table_Mass_0.2"},
            {"Table_Mass_0.3"},
            {"Table_Mass_0.4"},
            {"Table_Mass_0.5"},
            {"Table_Mass_0.6"},
            {"Table_Mass_0.7"},
            {"Table_Mass_0.8"},
            {"Table_Mass_0.9"},
            {"Table_Mass_1.0"},
            {"Table_Mass_1.1"},
            {"Table_Mass_1.2"},
            {"Table_Mass_1.3"}
        };
    }
};

#endif
