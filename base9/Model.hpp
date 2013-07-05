#ifndef MODEL_HPP
#define MODEL_HPP

#include <memory>

#include "MsRgbModel.hpp"
#include "Settings.hpp"
#include "ChabMsModel.hpp"
#include "DsedMsModel.hpp"
#include "GirardiMsModel.hpp"
#include "wdCooling.hpp"
#include "YaleMsModel.hpp"

using std::shared_ptr;

/*** Define a structure model that houses information about the evolution model ***/
class Model
{
  public:
    Model(shared_ptr<MsRgbModel> msRgbModel, shared_ptr<WdCoolingModel> wdCool, MsFilterSet filterSet)
        : mainSequenceEvol(msRgbModel), WDcooling(wdCool), filterSet(filterSet)
    {;}

    shared_ptr<MsRgbModel> mainSequenceEvol;
    shared_ptr<WdCoolingModel> WDcooling;

    MsFilterSet filterSet;

    int evoModel = 0;
    int brownDwarfEvol = 0;
    int IFMR = 0;
    int WDatm = 0;
    int numFilts = 0;
};

const Model makeModel(Settings &);

#endif
