#include <cstdio>

#include "Cluster.hpp"
#include "structures.hpp"

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
