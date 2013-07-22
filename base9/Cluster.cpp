#include <cstdio>

#include "Cluster.hpp"
#include "structures.hpp"

double Cluster::getAge() const
{
    return age;
}

double Cluster::getY() const
{
    return yyy;
}

double Cluster::getFeH() const
{
    return feh;
}

double Cluster::getMod() const
{
    return mod;
}

double Cluster::getAbs() const
{
    return abs;
}

double Cluster::getIfmrIntercept() const
{
    return ifmrIntercept;
}

double Cluster::getIfmrSlope() const
{
    return ifmrSlope;
}

double Cluster::getIfmrQuadCoef() const
{
    return ifmrQuadCoef;
}

void Cluster::setAge (const double newAge)
{
    age = newAge;
}

void Cluster::setY (const double newY)
{
    yyy = newY;
}

void Cluster::setFeH (const double newFeH)
{
    feh = newFeH;
}

void Cluster::setMod (const double newMod)
{
    mod = newMod;
}

void Cluster::setAbs (const double newAbs)
{
    abs = newAbs;
}

void Cluster::setIfmrIntercept (const double newIfmrIntercept)
{
    ifmrIntercept = newIfmrIntercept;
}

void Cluster::setIfmrSlope (const double newIfmrSlope)
{
    ifmrSlope = newIfmrSlope;
}

void Cluster::setIfmrQuadCoef (const double newIfmrQuadCoef)
{
    ifmrQuadCoef = newIfmrQuadCoef;
}

double Cluster::getParam(int p) const
{
    double temp;

    switch(p)
    {
        case AGE:
            temp = getAge();
            break;
        case YYY:
            temp = getY();
            break;
        case FEH:
            temp = getFeH();
            break;
        case MOD:
            temp = getMod();
            break;
        case ABS:
            temp = getAbs();
            break;
        case IFMR_INTERCEPT:
            temp = getIfmrIntercept();
            break;
        case IFMR_SLOPE:
            temp = getIfmrSlope();
            break;
        case IFMR_QUADCOEF:
            temp = getIfmrQuadCoef();
            break;
        default:
            throw std::out_of_range("Saving param matrix");
    }

    return temp;
}

void Cluster::setParam(int p, double v)
{
    switch(p)
    {
        case AGE:
            setAge(v);
            break;
        case YYY:
            setY(v);
            break;
        case FEH:
            setFeH(v);
            break;
        case MOD:
            setMod(v);
            break;
        case ABS:
            setAbs(v);
            break;
        case IFMR_INTERCEPT:
            setIfmrIntercept(v);
            break;
        case IFMR_SLOPE:
            setIfmrSlope(v);
            break;
        case IFMR_QUADCOEF:
            setIfmrQuadCoef(v);
            break;
        default:
            throw std::out_of_range("Saving param matrix");
    }
}
