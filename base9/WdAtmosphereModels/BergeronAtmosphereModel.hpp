#ifndef BERGATMOS_HPP
#define BERGATMOS_HPP

#include <array>
#include <map>
#include <string>
#include <utility>

#include <constants.

class BergeronAtmosphereModel : public WdCoolingModel
{
  public:
    virtual ~BergeronAtmosphereModel() {}

    virtual void loadModel (std::string path, MsFilter filterSet);
    virtual std::array<double, FILTS> teffToMags  (double wdLogTeff, double wdLogG, WdAtmosphere wdType);

  protected:
    std::map<double, std::vector<struct record>> hCurves;
    std::map<double, std::vector<struct record>> heCurves;
};
#endif
