#ifndef MODEL_HPP
#define MODEL_HPP

#include <memory>

#include "MsRgbModel.hpp"
#include "Settings.hpp"
#include "ChabMsModel.hpp"
#include "DsedMsModel.hpp"
#include "GirardiMsModel.hpp"
#include "YaleMsModel.hpp"

using std::shared_ptr;

/*** Define a structure model that houses information about the evolution model ***/
class Model
{
  public:
    Model(shared_ptr<MsRgbModel> msRgbModel, MsFilterSet filterSet)
        : mainSequenceEvol(msRgbModel), filterSet(filterSet)
    {;}

    int evoModel = 0;
    int brownDwarfEvol = 0;
    shared_ptr<MsRgbModel> mainSequenceEvol;
    MsFilterSet filterSet;
    int IFMR = 0;
    int WDcooling = 0;
    int WDatm = 0;
    int numFilts = 0;
};

const Model makeModel(Settings &);

#endif
