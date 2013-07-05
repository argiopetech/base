#include <iostream>
#include <memory>

#include "constants.hpp"
#include "Model.hpp"
#include "MsRgbModels/ChabMsModel.hpp"
#include "MsRgbModels/DsedMsModel.hpp"
#include "MsRgbModels/GirardiMsModel.hpp"
#include "MsRgbModels/YaleMsModel.hpp"
#include "WdCoolingModels/AlthausWdModel.hpp"
#include "WdCoolingModels/MontgomeryWdModel.hpp"
#include "WdCoolingModels/RenedoWdModel.hpp"
#include "WdCoolingModels/WoodWdModel.hpp"

using std::cerr;
using std::endl;
using std::shared_ptr;

const Model makeModel(Settings &s)
{
    shared_ptr<MsRgbModel> msModel;
    shared_ptr<WdCoolingModel> wdModel;
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

    switch (s.whiteDwarf.wdModel)
    {
        case WdModel::WOOD:
            wdModel = shared_ptr<WoodWdModel>(new WoodWdModel);
            break;
        case WdModel::MONTGOMERY:
            wdModel = shared_ptr<MontgomeryWdModel>(new MontgomeryWdModel);
            break;
        case WdModel::ALTHAUS:
            wdModel = shared_ptr<AlthausWdModel>(new AlthausWdModel);
            break;
        case WdModel::RENEDO:
            wdModel = shared_ptr<RenedoWdModel>(new RenedoWdModel);
            break;
        default:
            cerr << "***Error: No model found for white dwarf filter set " << static_cast<int>(s.whiteDwarf.wdModel) << ".***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (1);
    }

    Model model(msModel, wdModel, filterSet);

// !!! FIX ME !!!

    model.IFMR = s.whiteDwarf.ifmr;

    model.brownDwarfEvol = s.brownDwarf.bdModel;

    model.WDatm = BERGERON;

// END FIX ME

    return model;
}
