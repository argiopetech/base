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

    virtual void loadModel (std::string path, FilterSetName filterSet);
    virtual std::vector<double> teffToMags  (double wdLogTeff, double wdMass, WdAtmosphere wdType) const;

    virtual bool isSupported(FilterSetName) const;
};

#endif
