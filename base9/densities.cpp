#include <vector>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Cluster.hpp"
#include "Star.hpp"

#include "constants.hpp"
#include "densities.hpp"
#include "WhiteDwarf.hpp"

#include <iostream>

using std::vector;

constexpr double sqr(double a)
{
    return a * a;
}

static_assert(M_PI > 3.141592, "M_PI is defined and and at least 3.141592");
static_assert(M_PI < 3.15, "M_PI is defined and less than 3.15");

const double mf_sigma = 0.67729, mf_mu = -1.02;
const double loglog10 = log (log (10));

double Cluster::logPriorMass (double primaryMass) const
// Compute log prior density
{
    assert (primaryMass <= M_wd_up);

    double log_m1 = log10 (primaryMass);
    double logPrior = getLogMassNorm() + -0.5 * sqr (log_m1 - mf_mu) / sqr (mf_sigma) - log (primaryMass) - loglog10;

    return logPrior;
}

// Compute log prior density for cluster properties
double Cluster::logPrior (const Model &evoModels) const
{
    if ((age < evoModels.mainSequenceEvol->getMinAge())
        || (age > evoModels.mainSequenceEvol->getMaxAge())
        || (ifmrSlope < 0.0)
        || (abs < 0.0)
        || ((evoModels.IFMR == 11) && (ifmrQuadCoef < 0.0)))
    {
        throw InvalidCluster("Invalid cluster parameter in Cluster::logPrior");
    }

    // enforce monotonicity in IFMR
    if (evoModels.IFMR == 10)
    {
        double massLower = 0.15;
        double massUpper = M_wd_up;
        double massShift = 3.0;
        double angle = atan (ifmrSlope);
        double aa = cos (angle) * (1 + ifmrSlope * ifmrSlope);
        double xLower = aa * (massLower - massShift);
        double xUpper = aa * (massUpper - massShift);

        double dydx_xLower = ifmrQuadCoef * (xLower - xUpper);
        double dydx_xUpper = -dydx_xLower;

        double slopeLower = tan (angle + atan (dydx_xLower));
        double slopeUpper = tan (angle + atan (dydx_xUpper));

        // if IFMR is decreasing at either endpoint, reject
        if (slopeLower < 0.0 || slopeUpper < 0.0)
            throw std::range_error("IFMR decreasing at endpoint");
    }

    double prior = 0.0;
    //DS: with a uniform prior on carbonicity, the above won't change since log(1) = 0.

    if (priorVar[FEH] > EPSILON)
        prior += (-0.5) * sqr (feh - priorMean[FEH]) / priorVar[FEH];
    if (priorVar[MOD] > EPSILON)
        prior += (-0.5) * sqr (mod - priorMean[MOD]) / priorVar[MOD];
    if (priorVar[ABS] > EPSILON)
        prior += (-0.5) * sqr (abs - priorMean[ABS]) / priorVar[ABS];
    if (priorVar[YYY] > EPSILON)
        prior += (-0.5) * sqr (yyy - priorMean[YYY]) / priorVar[YYY];

    return prior;
}

static double scaledLogLike (const vector<double> &obsPhot, const vector<double> &variance, const vector<double> &mags, double varScale)
{
    double likelihood = 0.0;

    auto obsSize = obsPhot.size();

    for (size_t f = 0; f < obsSize; ++f)
    {
        auto varF = variance[f];

        if (varF > 1e-9)
        {
            varF *= varScale;
            likelihood -= 0.5 * (log (2 * M_PI * varF) + (sqr (mags[f] - obsPhot[f]) / varF));
        }
    }

    return likelihood;
}


// Calculates the posterior density for a stellar system
double StellarSystem::logPost (const Cluster &clust, const Model &evoModels, const Isochrone &isochrone)
{
    const vector<double> mags = deriveCombinedMags(clust, evoModels, isochrone);

    double logPrior = clust.logPriorMass (primary.mass);

    double likelihood = scaledLogLike (obsPhot, variance, mags, clust.varScale);

    return (logPrior + likelihood);
}
