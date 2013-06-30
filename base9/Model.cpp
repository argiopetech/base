#include <iostream>
#include <memory>

#include "constants.hpp"
#include "Model.hpp"
#include "MsRgbModel.hpp"

using std::cerr;
using std::endl;
using std::shared_ptr;

Model makeModel(Settings &s)
{
    shared_ptr<MsRgbModel> msModel;

    // !!! FIX ME !!!
    switch (s.mainSequence.msRgbModel) //    evoModels.mainSequenceEvol = settings.mainSequence.msRgbModel;
    {
        case MsModel::GIRARDI:
            msModel = shared_ptr<GirardiMsModel>(new GirardiMsModel);
            break;
        case MsModel::CHABHELIUM:
            msModel = shared_ptr<ChabMsModel>(new ChabMsModel);
            break;
        case MsModel::YALE:
            msModel = shared_ptr<YaleMsModel>(new YaleMsModel);
            break;
        case MsModel::DSED:
            msModel = shared_ptr<DsedMsModel>(new DsedMsModel);
            break;
        default:
            cerr << "***Error: No models found for main sequence evolution model " << static_cast<int>(s.mainSequence.msRgbModel) << ".***" << endl;
            cerr << "[Exiting...]\n" << endl;
            exit(1);
    }

    Model model(msModel);

    return model;
}
