#ifndef MODEL_HPP
#define MODEL_HPP

#include <memory>

#include "Filters.hpp"
#include "MsRgbModel.hpp"
#include "Settings.hpp"
#include "WdCoolingModel.hpp"
#include "WdAtmosphereModel.hpp"

using std::shared_ptr;

/*** Define a structure model that houses information about the evolution model ***/
class Model
{
  public:
    Model(shared_ptr<MsRgbModel> msRgbModel, shared_ptr<WdCoolingModel> wdCool, shared_ptr<WdAtmosphereModel> wdAtmos)
        : mainSequenceEvol(msRgbModel), WDcooling(wdCool), WDAtmosphere(wdAtmos)
    {;}

    void restrictFilters(const std::vector<std::string> &filters)
    {
        absCoeffs = Filters::calcAbsCoeffs(filters);
        mainSequenceEvol->restrictToFilters(filters);
        WDAtmosphere->restrictToFilters(filters);
    }

    shared_ptr<MsRgbModel> mainSequenceEvol;
    shared_ptr<WdCoolingModel> WDcooling;
    shared_ptr<WdAtmosphereModel> WDAtmosphere;

    std::vector<double> absCoeffs;

    int evoModel = 0;
    int IFMR = 0;
};

const Model makeModel(const Settings &);

#endif
