#ifndef BERGATMOS_HPP
#define BERGATMOS_HPP

#include <array>
#include <map>
#include <string>
#include <utility>

#include "constants.hpp"
#include "WdAtmosphereModel.hpp"

class BergeronAtmosphereModel : public WdAtmosphereModel
{
  public:
    BergeronAtmosphereModel()
    {
        availableFilters = {"U", "B", "V", "R", "I", "J", "H", "Ks", "u", "g", "r", "i", "z"
                           ,"y", "F435W", "F475W", "F555W", "F606W", "F625W", "F658N"
                           , "F775W", "F814W", "F850LP"};

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
            {"Table_Mass_1.2"}
        };
    }

    virtual ~BergeronAtmosphereModel() {}

    virtual void loadModel (std::string);

    virtual std::vector<double> teffToMags  (double, double, WdAtmosphere) const;
    virtual double teffToLogg (double, double, WdAtmosphere) const;

    virtual void restrictToFilters(const std::vector<std::string>&);

  protected:
    std::string dirName = "bergeron/";

    std::vector<std::string> files;
};

#endif
