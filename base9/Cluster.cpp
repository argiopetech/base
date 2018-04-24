#include <cmath>
#include <cstdio>

#include "Cluster.hpp"

double Cluster::getParam(int p) const
{
    double temp;

    switch(p)
    {
        case AGE:
            temp = age;
            break;
        case YYY:
            temp = yyy;
            break;
        case FEH:
            temp = feh;
            break;
        case MOD:
            temp = mod;
            break;
        case ABS:
            temp = abs;
            break;
        case CARBONICITY:
            temp = carbonicity;
            break;
        case IFMR_INTERCEPT:
            temp = ifmrIntercept;
            break;
        case IFMR_SLOPE:
            temp = ifmrSlope;
            break;
        case IFMR_QUADCOEF:
            temp = ifmrQuadCoef;
            break;
        default:
            throw std::out_of_range("Cluster::getParam()");
    }

    return temp;
}

void Cluster::setParam(int p, double v)
{
    switch(p)
    {
        case AGE:
            age = v;
            break;
        case YYY:
            yyy = v;
            break;
        case FEH:
            feh = v;
            break;
        case MOD:
            mod = v;
            break;
        case ABS:
            abs = v;
            break;
        case CARBONICITY:
            carbonicity = v;
            break;
        case IFMR_INTERCEPT:
            ifmrIntercept = v;
            break;
        case IFMR_SLOPE:
            ifmrSlope = v;
            break;
        case IFMR_QUADCOEF:
            ifmrQuadCoef = v;
            break;
        default:
            throw std::out_of_range("Cluster::setParam()");
    }
}

void Cluster::setM_wd_up(double M_wd_up)
{
    double p, q, c;
    double tup, tlow;

    constexpr double mf_sigma = 0.67729, mf_mu = -1.02;

    // computes normal distribution Phi(x) (integral from -Inf to x of normal density)
    // taken from: http://www.jstatsoft.org/v11/i04/v11i04.pdf
    auto Phi = [](double x)
    {
        long double s = x, t = 0, b = x, q = x * x, i = 1;

        while (s != t)
            s = (t = s) + (b *= q / (i += 2));

        return 0.5 + s * exp (-0.5 * q - 0.91893853320467274178L);
    };

    p = mf_mu + mf_sigma * mf_sigma * log (10);
    tup = (log10 (M_wd_up) - p) / (mf_sigma);
    tlow = (-1 - p) / mf_sigma;
    q = exp (-(mf_mu * mf_mu - p * p) / (2 * mf_sigma * mf_sigma));
    c = 1 / (q * mf_sigma * sqrt (2 * M_PI) * (Phi (tup) - Phi (tlow)));

    this->logMassNorm = log (c);
    this->M_wd_up = M_wd_up;
}
