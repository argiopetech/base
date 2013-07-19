#ifndef MODEL_HPP
#define MODEL_HPP

#include <memory>

#include "MsRgbModel.hpp"
#include "Settings.hpp"
#include "WdCoolingModel.hpp"
#include "MsFilterSet.hpp"

using std::shared_ptr;

/*** Define a structure model that houses information about the evolution model ***/
class Model
{
  public:
    Model(shared_ptr<MsRgbModel> msRgbModel, shared_ptr<MsFilterSet> filterSet, shared_ptr<WdCoolingModel> wdCool)
        : mainSequenceEvol(msRgbModel), filterSet(filterSet), WDcooling(wdCool)
    {;}

    shared_ptr<MsRgbModel> mainSequenceEvol;
    shared_ptr<MsFilterSet> filterSet;
    shared_ptr<WdCoolingModel> WDcooling;

    int evoModel = 0;
    int IFMR = 0;
    int WDatm = 0;
    int numFilts = 0;
};

const Model makeModel(const Settings &);

#endif
