#ifndef MODEL_HPP
#define MODEL_HPP

#include <memory>

#include "MsRgbModel.hpp"
#include "Settings.hpp"
#include "DsedMsModel.hpp"
#include "YaleMsModel.hpp"

using std::shared_ptr;

/*** Define a structure model that houses information about the evolution model ***/
class Model
{
  public:
    Model(shared_ptr<MsRgbModel> msRgbModel)
        : mainSequenceEvol(msRgbModel)
    {;}

    int evoModel = 0;
    int brownDwarfEvol = 0;
    shared_ptr<MsRgbModel> mainSequenceEvol;
    int IFMR = 0;
    int WDcooling = 0;
    int WDatm = 0;
    int filterSet = 0;
    int numFilts = 0;
};

Model makeModel(Settings &);

#endif
