#include <gtest/gtest.h>

#include "../LinearTransform.hpp"

#include <iostream>
#include <utility>

using std::make_pair;

TEST(LinearTransform, Interpolation)
{
    auto r1 = make_pair<double, double>(0, 10);
    auto r2 = make_pair<double, double>(0, 100);

    EXPECT_DOUBLE_EQ( 30.0, linearTransform<TransformMethod::Interp>(r1, r2,  3).val);
    EXPECT_DOUBLE_EQ(  0.0, linearTransform<TransformMethod::Interp>(r1, r2, -3).val);
    EXPECT_DOUBLE_EQ(100.0, linearTransform<TransformMethod::Interp>(r1, r2, 13).val);
}

TEST(LinearTransform, Extrapolation)
{
    auto r1 = make_pair<double, double>(0, 10);
    auto r2 = make_pair<double, double>(0, 100);

    EXPECT_DOUBLE_EQ( 30.0, linearTransform<TransformMethod::Extrap>(r1, r2,  3).val);
    EXPECT_DOUBLE_EQ(-30.0, linearTransform<TransformMethod::Extrap>(r1, r2, -3).val);
    EXPECT_DOUBLE_EQ(130.0, linearTransform<TransformMethod::Extrap>(r1, r2, 13).val);
}

TEST(LinearTransform, ExtrapLow)
{
    auto r1 = make_pair<double, double>(0, 10);
    auto r2 = make_pair<double, double>(0, 100);

    EXPECT_DOUBLE_EQ( 30.0, linearTransform<TransformMethod::ExtrapLow>(r1, r2,  3).val);
    EXPECT_DOUBLE_EQ(-30.0, linearTransform<TransformMethod::ExtrapLow>(r1, r2, -3).val);
    EXPECT_DOUBLE_EQ(100.0, linearTransform<TransformMethod::ExtrapLow>(r1, r2, 13).val);
}

TEST(LinearTransform, ExtrapHigh)
{
    auto r1 = make_pair<double, double>(0, 10);
    auto r2 = make_pair<double, double>(0, 100);

    EXPECT_DOUBLE_EQ( 30.0, linearTransform<TransformMethod::ExtrapHigh>(r1, r2,  3).val);
    EXPECT_DOUBLE_EQ(  0.0, linearTransform<TransformMethod::ExtrapHigh>(r1, r2, -3).val);
    EXPECT_DOUBLE_EQ(130.0, linearTransform<TransformMethod::ExtrapHigh>(r1, r2, 13).val);
}
