#include <gtest/gtest.h>

#include "../Model.hpp"
#include "../MsRgbModels/GirardiMsModel.hpp"
#include "../MsRgbModels/ChabMsModel.hpp"
#include "../MsRgbModels/YaleMsModel.hpp"
#include "../MsRgbModels/DsedMsModel.hpp"

#include <iostream>

TEST(Model, FilterSet)
{
    Settings s;
    s.whiteDwarf.wdModel = WdModel::WOOD;
    s.mainSequence.msRgbModel = MsModel::GIRARDI;

    s.mainSequence.filterSet = static_cast<FilterSetName>(5);
    EXPECT_EXIT(makeModel(s), ::testing::ExitedWithCode(1), "");
}

TEST(Model, makeModel)
{
    Settings s;
    s.whiteDwarf.wdModel = WdModel::WOOD;
    s.mainSequence.filterSet = FilterSetName::ACS;

    s.mainSequence.msRgbModel = MsModel::GIRARDI;
    EXPECT_EQ(typeid(*(new GirardiMsModel)), typeid(*makeModel(s).mainSequenceEvol));    

    s.mainSequence.msRgbModel = MsModel::CHABHELIUM;
    EXPECT_EQ(typeid(*(new ChabMsModel)), typeid(*makeModel(s).mainSequenceEvol));    

    s.mainSequence.msRgbModel = MsModel::YALE;
    EXPECT_EQ(typeid(*(new YaleMsModel)), typeid(*makeModel(s).mainSequenceEvol));

    s.mainSequence.msRgbModel = MsModel::DSED;
    EXPECT_EQ(typeid(*(new DsedMsModel)), typeid(*makeModel(s).mainSequenceEvol));
}
