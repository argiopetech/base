#include <gtest/gtest.h>

#include "../Model.hpp"

#include <iostream>

int verbose, needMassNow;
int *useFilt;
double *ltau, *aFilt;

TEST(Model, makeModel)
{
    Settings s;

    s.mainSequence.msRgbModel = MsModel::GIRARDI;
    EXPECT_EQ(typeid(*(new GirardiMsModel)), typeid(*makeModel(s).mainSequenceEvol));    

    s.mainSequence.msRgbModel = MsModel::CHABHELIUM;
    EXPECT_EQ(typeid(*(new ChabMsModel)), typeid(*makeModel(s).mainSequenceEvol));    

    s.mainSequence.msRgbModel = MsModel::YALE;
    EXPECT_EQ(typeid(*(new YaleMsModel)), typeid(*makeModel(s).mainSequenceEvol));

    s.mainSequence.msRgbModel = MsModel::DSED;
    EXPECT_EQ(typeid(*(new DsedMsModel)), typeid(*makeModel(s).mainSequenceEvol));
}
