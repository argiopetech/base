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
    MsFilterSet filterSet;

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

    switch (s.mainSequence.filterSet)
    {
        case MsFilterSet::UBVRIJHK:
        case MsFilterSet::ACS:
        case MsFilterSet::SDSS:
            filterSet = s.mainSequence.filterSet;
            break;
        default:
            cerr << "***Error: No models found for filter set " << static_cast<int>(s.mainSequence.filterSet) << ".***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (1);
    }

    Model model(msModel, filterSet);

    return model;
}
