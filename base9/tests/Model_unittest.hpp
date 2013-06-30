#include <gtest/gtest.h>

#include "../Model.hpp"
#include "../YaleMsModel.hpp"

#include <iostream>

int verbose;
int *useFilt;
double *ltau, *aFilt;

TEST(Model, makeModel)
{
    Settings s;

    

    s.mainSequence.msRgbModel = YALE;
    EXPECT_EQ(typeid(*(new YaleMsModel)), typeid(*makeModel(s).mainSequenceEvol));

    s.mainSequence.msRgbModel = DSED;
    EXPECT_EQ(typeid(*(new DsedMsModel)), typeid(*makeModel(s).mainSequenceEvol));
}
