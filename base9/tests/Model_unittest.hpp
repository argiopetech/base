#include <gtest/gtest.h>

#include "../Model.hpp"
#include "../MsRgbModels/GirardiMsModel.hpp"
#include "../MsRgbModels/ChabMsModel.hpp"
#include "../MsRgbModels/YaleMsModel.hpp"
#include "../MsRgbModels/DsedMsModel.hpp"

#include <iostream>

int *useFilt;
double *ltau, *aFilt;

TEST(Model, FilterSet)
{
    Settings s;
    s.whiteDwarf.wdModel = WOOD;
    s.mainSequence.msRgbModel = MsModel::GIRARDI;

    s.mainSequence.filterSet = MsFilterSet::UBVRIJHK;
    EXPECT_EQ(MsFilterSet::UBVRIJHK, makeModel(s).filterSet);

    s.mainSequence.filterSet = MsFilterSet::ACS;
    EXPECT_EQ(MsFilterSet::ACS, makeModel(s).filterSet);

    s.mainSequence.filterSet = MsFilterSet::SDSS;
    EXPECT_EQ(MsFilterSet::SDSS, makeModel(s).filterSet);

    s.mainSequence.filterSet = static_cast<MsFilterSet>(5);
    EXPECT_EXIT(makeModel(s), ::testing::ExitedWithCode(1), "");
}

TEST(Model, makeModel)
{
    Settings s;
    s.mainSequence.filterSet = MsFilterSet::ACS;

    s.mainSequence.msRgbModel = MsModel::GIRARDI;
    EXPECT_EQ(typeid(*(new GirardiMsModel)), typeid(*makeModel(s).mainSequenceEvol));    

    s.mainSequence.msRgbModel = MsModel::CHABHELIUM;
    EXPECT_EQ(typeid(*(new ChabMsModel)), typeid(*makeModel(s).mainSequenceEvol));    

    s.mainSequence.msRgbModel = MsModel::YALE;
    EXPECT_EQ(typeid(*(new YaleMsModel)), typeid(*makeModel(s).mainSequenceEvol));

    s.mainSequence.msRgbModel = MsModel::DSED;
    EXPECT_EQ(typeid(*(new DsedMsModel)), typeid(*makeModel(s).mainSequenceEvol));
}
