#include <iostream>
#include <memory>

#include "gBergMag.hpp"

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

using std::cout;
using std::cerr;
using std::endl;
using std::shared_ptr;

namespace internal
{
    shared_ptr<MsRgbModel> createMsRgbModel(MsModel model)
    {
        switch (model)
        {
            case MsModel::GIRARDI:
                return shared_ptr<GirardiMsModel>(new GirardiMsModel);
            case MsModel::CHABHELIUM:
                return shared_ptr<ChabMsModel>(new ChabMsModel);
            case MsModel::YALE:
                return shared_ptr<YaleMsModel>(new YaleMsModel);
            case MsModel::DSED:
                return shared_ptr<DsedMsModel>(new DsedMsModel);
            default:
                std::cerr << "***Error: No models found for main sequence evolution model " << static_cast<int>(model) << ".***" << std::endl;
                std::cerr << "[Exiting...]\n" << std::endl;
                exit(1);
        }
    }

    shared_ptr<MsFilterSet> createMsFilterSet(MsFilter filter)
    {
        switch (filter)
        {
            case MsFilter::UBVRIJHK:
                return shared_ptr<UBVRIJHK>(new UBVRIJHK);
            case MsFilter::ACS:
                return shared_ptr<ACS>(new ACS);
            case MsFilter::SDSS:
                return shared_ptr<SDSS>(new SDSS);
            default:
                cerr << "***Error: No models found for filter set " << static_cast<int>(filter) << ".***" << endl;
                cerr << "[Exiting...]" << endl;
                exit (1);
        }
    }

    shared_ptr<WdCoolingModel> createWdCoolingModel(WdModel model)
    {
        switch (model)
        {
            case WdModel::WOOD:
                return shared_ptr<WoodWdModel>(new WoodWdModel);
            case WdModel::MONTGOMERY:
                return shared_ptr<MontgomeryWdModel>(new MontgomeryWdModel);
            case WdModel::ALTHAUS:
                return shared_ptr<AlthausWdModel>(new AlthausWdModel);
            case WdModel::RENEDO:
                return shared_ptr<RenedoWdModel>(new RenedoWdModel);
            default:
                cerr << "***Error: No model found for white dwarf filter set " << static_cast<int>(model) << ".***" << endl;
                cerr << "[Exiting...]" << endl;
                exit (1);
        }
    }
}

const Model makeModel(const Settings &s)
{
    cout << "\nReading models..." << std::flush;

    Model model( internal::createMsRgbModel(s.mainSequence.msRgbModel)
               , internal::createMsFilterSet(s.mainSequence.filterSet)
               , internal::createWdCoolingModel(s.whiteDwarf.wdModel));

// !!! FIX ME !!!

    model.IFMR = s.whiteDwarf.ifmr;

    model.WDatm = BERGERON;

// END FIX ME

    model.mainSequenceEvol->loadModel(s.files.models, s.mainSequence.filterSet);
    model.WDcooling->loadModel(s.files.models);

    loadBergeron (s.files.models, s.mainSequence.filterSet);

    cout << " Done." << endl;

    return model;
}
