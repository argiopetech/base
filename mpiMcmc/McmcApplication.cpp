#include <iostream>

#include <random>

#include <cmath>

#include "McmcApplication.hpp"

using std::cerr;
using std::endl;

// Decides whether to accept a proposed cluster property
bool McmcApplication::acceptClustMarg (const double logPostCurr, const double logPostProp)
{
    if (isinf (logPostProp))
    {
//        cerr << "-Inf posterior proposed and rejected" << endl;
        rejected += 1;
        return false;
    }

    double alpha = logPostProp - logPostCurr;

    if (alpha >= 0)             // Short circuit exit to the MH algorithm
    {
        accepted += 1;
        return true;
    }

    double u = std::generate_canonical<double, 53>(gen);

    if (u < 1.e-15)
        u = 1.e-15;
    u = log (u);

    if (u < alpha)
    {
        accepted += 1;
        return true;
    }
    else
    {
        rejected += 1;
        return false;
    }
}
