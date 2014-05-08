#ifndef INVALIDATMOS_HPP
#define INVALIDATMOS_HPP

#include <array>
#include <map>
#include <string>
#include <utility>

#include "constants.hpp"
#include "WdAtmosphereModel.hpp"

class InvalidAtmosphereModel : public WdAtmosphereModel, public InvalidModel
{
  public:
    virtual ~InvalidAtmosphereModel() {}

    virtual std::vector<double> teffToMags  (double, double, WdAtmosphere) const
    {
        throw InvalidModelError("Called teffToMags() in invalid atmosphere model");
    }
};

#endif
