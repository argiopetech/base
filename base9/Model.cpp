#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "constants.hpp"
#include "Model.hpp"
#include "MsRgbModels/ChabMsModel.hpp"
#include "MsRgbModels/OldDsedMsModel.hpp"
#include "MsRgbModels/NewDsedMsModel.hpp"
#include "MsRgbModels/GirardiMsModel.hpp"
#include "MsRgbModels/InvalidMsModel.hpp"
#include "WdCoolingModels/AlthausWdModel.hpp"
#include "WdCoolingModels/MontgomeryWdModel.hpp"
#include "WdCoolingModels/RenedoWdModel.hpp"
#include "WdCoolingModels/WoodWdModel.hpp"
#include "WdAtmosphereModels/BergeronAtmosphereModel.hpp"
#include "WdAtmosphereModels/InvalidAtmosphereModel.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::shared_ptr;
using std::string;
using std::vector;

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
                cerr << "Yale models are not currently supported." << endl;
                exit(0);
            case MsModel::OLD_DSED:
                return shared_ptr<OldDsedMsModel>(new OldDsedMsModel);
            case MsModel::NEW_DSED:
                return shared_ptr<NewDsedMsModel>(new NewDsedMsModel);
            default:
                std::cerr << "***Error: No models found for main sequence evolution model " << static_cast<int>(model) << ".***" << std::endl;
                std::cerr << "[Exiting...]\n" << std::endl;
                exit(1);
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

    shared_ptr<WdAtmosphereModel> createWdAtmosphereModel(WdAtmosphereModelSet model)
    {
        switch (model)
        {
            case WdAtmosphereModelSet::BERGERON:
                return shared_ptr<BergeronAtmosphereModel>(new BergeronAtmosphereModel);
            default:
                cerr << "***Error: No model found for white dwarf atmosphere set " << static_cast<int>(model) << ".***" << endl;
                cerr << "[Exiting...]" << endl;
                exit (1);
        }
    }
}

const Model makeModel(const Settings &s)
{
    cout << "Reading models..." << std::flush;

    Model model( internal::createMsRgbModel(s.mainSequence.msRgbModel)
               , internal::createWdCoolingModel(s.whiteDwarf.wdModel)
               , internal::createWdAtmosphereModel(WdAtmosphereModelSet::BERGERON));
                 
    if (! model.mainSequenceEvol->isSupported(s.mainSequence.filterSet))
    {
        cout << "\n\tMSRGB model does not support the selected filter set. Only WD operations will be supported." << endl;
        model.mainSequenceEvol = shared_ptr<InvalidMsModel>(new InvalidMsModel());
    }

    if (! model.WDAtmosphere->isSupported(s.mainSequence.filterSet))
    {
        cout << "\n\tWD Atmosphere model does not support the selected filter set. Only MS operations will be supported." << endl;
        model.WDAtmosphere = shared_ptr<InvalidAtmosphereModel>(new InvalidAtmosphereModel());
    }

// !!! FIX ME !!!

    model.IFMR = s.whiteDwarf.ifmr;

// END FIX ME

    model.mainSequenceEvol->loadModel(s.files.models, s.mainSequence.filterSet);
    model.WDcooling->loadModel(s.files.models, s.mainSequence.filterSet);
    model.WDAtmosphere->loadModel(s.files.models, s.mainSequence.filterSet);

    cout << " Done.\n" << endl;

    return model;
}
