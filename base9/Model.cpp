#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "constants.hpp"
#include "Model.hpp"
#include "MsRgbModels/ChabMsModel.hpp"
#include "MsRgbModels/OldDsedMsModel.hpp"
#include "MsRgbModels/NewDsedMsModel.hpp"
#include "MsRgbModels/ParsecMsModel.hpp"
#include "MsRgbModels/GirardiMsModel.hpp"
#include "WdCoolingModels/AlthausWdModel.hpp"
#include "WdCoolingModels/MontgomeryWdModel.hpp"
#include "WdCoolingModels/RenedoWdModel.hpp"
#include "WdCoolingModels/WoodWdModel.hpp"
#include "WdAtmosphereModels/BergeronAtmosphereModel.hpp"

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
            case MsModel::PARSEC:
                return shared_ptr<ParsecMsModel>(new ParsecMsModel);
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
    if (s.verbose)
        cout << "Reading models..." << std::flush;

    Model model( internal::createMsRgbModel(s.mainSequence.msRgbModel)
               , internal::createWdCoolingModel(s.whiteDwarf.wdModel)
               , internal::createWdAtmosphereModel(WdAtmosphereModelSet::BERGERON)
               , s.verbose);

// !!! FIX ME !!!

    model.IFMR = s.whiteDwarf.ifmr;

// END FIX ME

    model.mainSequenceEvol->loadModel(s.files.models);
    model.WDcooling->loadModel(s.files.models);
    model.WDAtmosphere->loadModel(s.files.models);

    if (s.verbose)
        cout << " Done.\n" << endl;

    return model;
}
