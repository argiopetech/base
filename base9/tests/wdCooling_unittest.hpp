#include <gtest/gtest.h>

#include "../wdCooling.hpp"

#include <iostream>

TEST(WdCooling, wdMassToTeffAndRadius)
{
    double vars[6][5] = {
        {8.796, 8.6418, 0.711294, 8.90714, 4.23269},
        {8.796, 8.62552, 0.716454, 8.90406, 4.22367},
        {8.796, 8.60251, 0.723936, 8.89966, 4.21289},
        {8.796, 8.58781, 0.728838, 8.89681, 4.20705},
        {8.796, 8.42118, 0.79179, 8.86338, 4.17788},
        {8.796, 8.42424, 0.7905, 8.86406, 4.17788}};


    loadWDCool("/home/elliot/Projects/stellar_evolution/test/hyades2/models/", ALTHAUS);

    double logRadius = 0.0;

    for (int i = 0; i < 6; ++i)
    {
        ASSERT_NEAR(vars[i][4], wdMassToTeffAndRadius(vars[i][0], 0.36, vars[i][1], vars[i][2], &logRadius), 0.0001);
        ASSERT_NEAR(vars[i][3], logRadius, 0.0001);
    }
}
