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
        availableFilters = {"U", "B", "V", "R", "I", "J", "H", "K", "u", "g", "r", "i", "z"};
    }

    virtual ~BergeronAtmosphereModel() {}

    virtual void loadModel (std::string);
    virtual std::vector<double> teffToMags  (double, double, WdAtmosphere) const;
    virtual void restrictToFilters(const std::vector<std::string>&);
};

#endif
