#include <array>
#include <vector>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Cluster.hpp"
#include "Star.hpp"

#include "constants.hpp"
#include "structures.hpp"
#include "densities.hpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::vector;

constexpr double sqr(double a)
{
    return a * a;
}

static_assert(M_PI > 3.141592, "M_PI is defined and and at least 3.141592");
static_assert(M_PI < 3.15, "M_PI is defined and less than 3.15");

const double mf_sigma = 0.67729, mf_mu = -1.02;
const double loglog10 = log (log (10));

// Forward Declarations
double Phi (double);

namespace evil
{
    // Calculate the mass normalization factor once so that we don't have to
    // calculate it every time we want the mass prior
    class logMassNorm
    {
      private:
        logMassNorm(double M_wd_up)
        {
            double p, q, c;
            double tup, tlow;

            p = mf_mu + mf_sigma * mf_sigma * log (10);
            tup = (log10 (M_wd_up) - p) / (mf_sigma);
            tlow = (-1 - p) / mf_sigma;
            q = exp (-(mf_mu * mf_mu - p * p) / (2 * mf_sigma * mf_sigma));
            c = 1 / (q * mf_sigma * sqrt (2 * M_PI) * (Phi (tup) - Phi (tlow)));

            this->var = log (c);
        }
        logMassNorm(const logMassNorm&) = delete;
        logMassNorm(const logMassNorm&&) = delete;
        logMassNorm& operator=(const logMassNorm&) = delete;

        double var;

      public:
        static logMassNorm& getInstance(const Cluster &clust)
        {
            static logMassNorm instance(clust.M_wd_up);

            return instance;
        }

        double getLogMassNorm() const { return var; }
    };
}

double logPriorMass (const StellarSystem &system, const Cluster &clust)
// Compute log prior density
{
    double mass1, log_m1, logPrior = 0.0;

    mass1 = system.primary.mass;

    assert (mass1 > 0.1 && mass1 <= clust.M_wd_up);

    log_m1 = log10 (mass1);
    logPrior = evil::logMassNorm::getInstance(clust).getLogMassNorm() + -0.5 * sqr (log_m1 - mf_mu) / (sqr (mf_sigma)) - log (mass1) - loglog10;
    return logPrior;
}

// Compute log prior density for cluster properties
double logPriorClust (const Cluster &clust, const Model &evoModels)
{
    if ((clust.age < evoModels.mainSequenceEvol->getMinAge())
        || (clust.age > evoModels.mainSequenceEvol->getMaxAge())
        || (clust.ifmrSlope < 0.0)
        || (clust.abs < 0.0)
        || ((evoModels.IFMR == 11) && (clust.ifmrQuadCoef < 0.0)))
    {
        throw InvalidCluster("Invalid cluster parameter");
    }

    // enforce monotonicity in IFMR
    if (evoModels.IFMR == 10)
    {
        double massLower = 0.15;
        double massUpper = clust.M_wd_up;
        double massShift = 3.0;
        double angle = atan (clust.ifmrSlope);
        double aa = cos (angle) * (1 + clust.ifmrSlope * clust.ifmrSlope);
        double xLower = aa * (massLower - massShift);
        double xUpper = aa * (massUpper - massShift);

        double dydx_xLower = clust.ifmrQuadCoef * (xLower - xUpper);
        double dydx_xUpper = -dydx_xLower;

        double slopeLower = tan (angle + atan (dydx_xLower));
        double slopeUpper = tan (angle + atan (dydx_xUpper));

        // if IFMR is decreasing at either endpoint, reject
        if (slopeLower < 0.0 || slopeUpper < 0.0)
            throw std::range_error("IFMR decreasing at endpoint");
    }

    double prior = 0.0;
    //DS: with a uniform prior on carbonicity, the above won't change since log(1) = 0.

    if (clust.priorVar.at(FEH) > EPSILON)
        prior += (-0.5) * sqr (clust.feh - clust.priorMean.at(FEH)) / clust.priorVar.at(FEH);
    if (clust.priorVar.at(MOD) > EPSILON)
        prior += (-0.5) * sqr (clust.mod - clust.priorMean.at(MOD)) / clust.priorVar.at(MOD);
    if (clust.priorVar.at(ABS) > EPSILON)
        prior += (-0.5) * sqr (clust.abs - clust.priorMean.at(ABS)) / clust.priorVar.at(ABS);
    if (clust.priorVar.at(YYY) > EPSILON)
        prior += (-0.5) * sqr (clust.yyy - clust.priorMean.at(YYY)) / clust.priorVar.at(YYY);

    return prior;
}

double logLikelihood (int numFilts, const StellarSystem &system, const array<double, FILTS> &mags)
// Computes log likelihood
{
    int i;
    double likelihood = 0.0;

    for (i = 0; i < numFilts; i++)
    {
        if (system.variance.at(i) > 1e-9)
            likelihood -= 0.5 * (log (2 * M_PI * system.variance.at(i)) + (sqr (mags.at(i) - system.obsPhot.at(i)) / system.variance.at(i)));
    }

    return likelihood;
}

double tLogLikelihood (int numFilts, const StellarSystem &system, const array<double, FILTS> &mags)
// Computes log likelihood
{
    int i;
    double likelihood = 0.0;
    double dof = 3.0;
    double quadratic_sum = 0.0;

    for (i = 0; i < numFilts; i++)
    {
        if (system.variance.at(i) > 1e-9)
        {
            quadratic_sum += sqr (mags.at(i) - system.obsPhot.at(i)) / system.variance.at(i);
            likelihood -= 0.5 * (log (M_PI * system.variance.at(i)));
        }
    }
    likelihood += lgamma (0.5 * (dof + (double) numFilts)) - lgamma (0.5 * dof);
    likelihood -= 0.5 * (dof + (double) numFilts) * log (1 + quadratic_sum);

    return likelihood;
}

double scaledLogLike (int numFilts, const StellarSystem &system, double varScale, const array<double, FILTS> &mags)
// Computes log likelihood
{
    int i;
    double likelihood = 0.0;

    for (i = 0; i < numFilts; i++)
    {
        if (system.variance.at(i) > 1e-9)
        {
            likelihood -= 0.5 * (log (2 * M_PI * varScale * system.variance.at(i)) + (sqr (mags.at(i) - system.obsPhot.at(i)) / (varScale * system.variance.at(i))));
        }

    }
    return likelihood;
}


double logPost1Star (const StellarSystem &system, const Cluster &clust, const Model &evoModels, const vector<int> &filters)
// Compute posterior density for 1 star:
{
    // AGBt_zmass never set because age and/or metallicity out of range of models.
    if (clust.AGBt_zmass < EPS)
    {
        throw WDBoundsError("Bounds error in evolve.cpp");
    }

    double likelihood = 0.0, logPrior = 0.0;

    const array<double, FILTS> mags = system.deriveCombinedMags(clust, evoModels, filters);

    logPrior = logPriorMass (system, clust);

    likelihood = scaledLogLike (evoModels.numFilts, system, clust.varScale, mags);

    return (logPrior + likelihood);
}


// computes normal distribution Phi(x) (integral from -Inf to x of normal density)
// taken from: http://www.jstatsoft.org/v11/i04/v11i04.pdf
double Phi (double x)
{
    long double s = x, t = 0, b = x, q = x * x, i = 1;

    while (s != t)
        s = (t = s) + (b *= q / (i += 2));

    return 0.5 + s * exp (-0.5 * q - 0.91893853320467274178L);
}

// do not call this routine with nu = 2, which wouldn't make much sense anyway
double logTDens (double x, double mean, double var, double nu)
{
    double logp = 0;
    double s;

    s = sqrt (nu / (var * (nu - 2)));

    logp = log (s) + gamma(nu) - 3.5 * log (1 + pow (s * (x - mean), 2) / nu);

    return logp;
}
