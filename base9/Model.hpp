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

using std::cerr;
using std::cout;
using std::endl;
using std::shared_ptr;

/*** Define a structure model that houses information about the evolution model ***/
class Model
{
  public:
    Model(shared_ptr<MsRgbModel> msRgbModel, shared_ptr<WdCoolingModel> wdCool, shared_ptr<WdAtmosphereModel> wdAtmos, bool verbose)
        : mainSequenceEvol(msRgbModel), WDcooling(wdCool), WDAtmosphere(wdAtmos), verbose(verbose)
    {;}

    void restrictFilters(const std::vector<std::string> &filters, bool allowInvalid)
    {
        absCoeffs = Filters::calcAbsCoeffs(filters);

        try
        {
            mainSequenceEvol->restrictToFilters(filters, allowInvalid);
        }
        catch (InvalidModelError &e)
        {
            if (verbose)
            {
                cerr << "\n***Warning: \"" << e.what() << "\" is not available in the selected MSRGB model\n"
                     << "            This is non-fatal if you aren't using the MSRGB models" << endl;
            }
            mainSequenceEvol.reset(new InvalidMsModel());
        }

        try
        {
            WDAtmosphere->restrictToFilters(filters);
        }
        catch (InvalidModelError &e)
        {
            if (verbose)
            {
                cerr << "\n***Warning: \"" << e.what() << "\" is not available in the selected WD Atmosphere model\n"
                     << "            This is non-fatal if you aren't using the WD models" << endl;
            }

            WDAtmosphere.reset(new InvalidAtmosphereModel());
        }
    }

    shared_ptr<MsRgbModel> mainSequenceEvol;
    shared_ptr<WdCoolingModel> WDcooling;
    shared_ptr<WdAtmosphereModel> WDAtmosphere;

    std::vector<double> absCoeffs;

    int evoModel = 0;
    int IFMR = 0;

  private:
    bool verbose;
};

const Model makeModel(const Settings &);

#endif
