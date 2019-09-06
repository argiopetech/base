#ifndef BERG2019ATMOS_HPP
#define BERG2019ATMOS_HPP

#include <array>
#include <map>
#include <string>
#include <utility>

#include "constants.hpp"
#include "WdAtmosphereModel.hpp"

class Bergeron2019AtmosphereModel : public BergeronAtmosphereModel
{
  public:
    Bergeron2019AtmosphereModel()
    {
        dirName = "bergeron_2019/";

        availableFilters = { "U", "B", "V", "R", "I",                // Johnson-Kron-Cousins UBVRI
                             "J", "H", "Ks",                         // 2MASS JHKs, where Ks means "K short"
                             "Y_MKO", "J_MKO", "H_MKO", "K_MKO",     // Mauna Kea Observatory (MKO) YJHK
                             "W1", "W2", "W3", "W4",                 // Wide-field Infrared Survey Explorer (WISE) W1 to W4
                             "S3.6", "S4.5", "S5.8", "S8.0",         // Spitzer IRAC 3.6, 4.5, 5.8, 8.0 microns
                             "u", "g", "r", "i", "z",                // SDSS
                             "g_ps", "r_ps", "i_ps", "z_ps", "y_ps", // Pan-STARRS
                             "G", "G_BP", "G_RP",                    // Gaia
                             "FUV", "NUV" };                         // Galex
    }
};

#endif
