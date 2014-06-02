#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <memory>

#include "Filters.hpp"
#include "MsRgbModel.hpp"
#include "MsRgbModels/InvalidMsModel.hpp"
#include "Settings.hpp"
#include "WdCoolingModel.hpp"
#include "WdAtmosphereModel.hpp"
#include "WdAtmosphereModels/InvalidAtmosphereModel.hpp"

using std::cout;
using std::endl;
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

        try
        {
            mainSequenceEvol->restrictToFilters(filters);
        }
        catch (InvalidModelError &e)
        {
            cout << "Couldn't find filter \"" << e.what() << "\" in selected MSRGB model" << endl;
            mainSequenceEvol.reset(new InvalidMsModel());
        }

        try
        {
            WDAtmosphere->restrictToFilters(filters);
        }
        catch (InvalidModelError &e)
        {
            cout << "Couldn't find filter \"" << e.what() << "\" in selected WD Atmosphere model" << endl;
            WDAtmosphere.reset(new InvalidAtmosphereModel());
        }
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
