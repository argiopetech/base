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
    }

    virtual ~BergeronAtmosphereModel() {}

    virtual void loadModel (std::string);

    virtual std::vector<double> teffToMags  (double, double, WdAtmosphere) const;
    virtual double teffToLogg (double, double, WdAtmosphere) const;

    virtual void restrictToFilters(const std::vector<std::string>&);

  protected:
    std::string dirName = "bergeron/";
};

#endif
