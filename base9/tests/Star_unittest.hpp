#include <gtest/gtest.h>

#include "../Matrix.hpp"
#include "../Star.hpp"

TEST(Star, Constructor)
{
    Star s;

    // Ensures the matrix fill worked
    // Failure case: (auto a) instead of (auto &a)
    EXPECT_TRUE([&] {
            bool b = true;

            for (auto a : s.beta)
                for (auto e : a)
                    b &= (e == 0.0);

            return b;
        }());
}
