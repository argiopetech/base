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
    virtual ~BergeronAtmosphereModel() {}

    virtual void loadModel (std::string path, MsFilter filterSet);
    virtual std::array<double, FILTS> teffToMags  (double wdLogTeff, double wdLogG, WdAtmosphere wdType) const;
};
#endif
