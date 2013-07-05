#include <gtest/gtest.h>

#include "../wdCooling.hpp"

#include <iostream>

TEST(WdCooling, wdMassToTeffAndRadius)
{
    loadWDCool("/home/elliot/Projects/stellar_evolution/test/hyades2/models/", ALTHAUS);

    double logRadius = 0.0;

    ASSERT_NEAR(4.23269, wdMassToTeffAndRadius(8.796, 0.36, 8.6418, 0.711294, logRadius), 0.0001);
    ASSERT_NEAR(8.90714, logRadius, 0.0001);

    ASSERT_NEAR(4.22367, wdMassToTeffAndRadius(8.796, 0.36, 8.62552, 0.716454, logRadius), 0.0001);
    ASSERT_NEAR(8.90406, logRadius, 0.0001);

    ASSERT_NEAR(4.21289, wdMassToTeffAndRadius(8.796, 0.36, 8.60251, 0.723936, logRadius), 0.0001);
    ASSERT_NEAR(8.89966, logRadius, 0.0001);

    ASSERT_NEAR(4.20706, wdMassToTeffAndRadius(8.796, 0.36, 8.58781, 0.728838, logRadius), 0.0001);
    ASSERT_NEAR(8.89681, logRadius, 0.0001);

    ASSERT_NEAR(4.17788, wdMassToTeffAndRadius(8.796, 0.36, 8.42118, 0.79179, logRadius), 0.0001);
    ASSERT_NEAR(8.86338, logRadius, 0.0001);

    ASSERT_NEAR(4.17788, wdMassToTeffAndRadius(8.796, 0.36, 8.42424, 0.7905, logRadius), 0.0001);
    ASSERT_NEAR(8.86406, logRadius, 0.0001);
}
