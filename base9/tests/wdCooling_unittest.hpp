#include <gtest/gtest.h>

#include "../WdCoolingModels/AlthausWdModel.hpp"

#include <iostream>
#include <utility>

using std::make_pair;

TEST(WdCooling, wdMassToTeffAndRadius)
{
    AlthausWdModel wdCool;

    wdCool.loadModel("/home/elliot/Projects/stellar_evolution/test/hyades2/models/");

    ASSERT_NEAR((make_pair<double, double>(4.23269, 8.90714)).first, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.6418, 0.711294).first, 0.0001);
    ASSERT_NEAR((make_pair<double, double>(4.23269, 8.90714)).second, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.6418, 0.711294).second, 0.0001);

    ASSERT_NEAR((make_pair<double, double>(4.22367, 8.90406)).first, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.62552, 0.716454).first, 0.0001);
    ASSERT_NEAR((make_pair<double, double>(4.22367, 8.90406)).second, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.62552, 0.716454).second, 0.0001);

    ASSERT_NEAR((make_pair<double, double>(4.21289, 8.89966)).first, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.60251, 0.723936).first, 0.0001);
    ASSERT_NEAR((make_pair<double, double>(4.21289, 8.89966)).second, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.60251, 0.723936).second, 0.0001);

    ASSERT_NEAR((make_pair<double, double>(4.20706, 8.89681)).first, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.58781, 0.728838).first, 0.0001);
    ASSERT_NEAR((make_pair<double, double>(4.20706, 8.89681)).second, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.58781, 0.728838).second, 0.0001);

    ASSERT_NEAR((make_pair<double, double>(4.17788, 8.86338)).first, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.42118, 0.79179).first, 0.0001);
    ASSERT_NEAR((make_pair<double, double>(4.17788, 8.86338)).second, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.42118, 0.79179).second, 0.0001);

    ASSERT_NEAR((make_pair<double, double>(4.17788, 8.86406)).first, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.42424, 0.7905).first, 0.0001);
    ASSERT_NEAR((make_pair<double, double>(4.17788, 8.86406)).second, wdCool.wdMassToTeffAndRadius(8.796, 0.36, 8.42424, 0.7905).second, 0.0001);
}
