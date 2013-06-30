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
        // case GIR:
        //      break;
        // case CHABHELIUM:
        //      break;
        case YALE:
            msModel = shared_ptr<YaleMsModel>(new YaleMsModel);
            break;
        case DSED:
            msModel = shared_ptr<DsedMsModel>(new DsedMsModel);
            break;
        default:
            cerr << "***Error: No models found for main sequence evolution model " << s.mainSequence.msRgbModel << ".***" << endl;
            cerr << "[Exiting...]\n" << endl;
            exit(1);
    }

    Model model(msModel);

    return model;
}
