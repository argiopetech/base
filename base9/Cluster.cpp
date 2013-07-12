#include <cstdio>

#include "Cluster.hpp"
#include "structures.hpp"

double Cluster::getAge() const
{
    return parameter[AGE];
}

double Cluster::getY() const
{
    return parameter[YYY];
}

double Cluster::getFeH() const
{
    return parameter[FEH];
}

double Cluster::getMod() const
{
    return parameter[MOD];
}

double Cluster::getAbs() const
{
    return parameter[ABS];
}

void Cluster::setAge (const double newAge)
{
    parameter[AGE] = newAge;
}

void Cluster::setY (const double newY)
{
    parameter[YYY] = newY;
}

void Cluster::setFeH (const double newFeH)
{
    parameter[FEH] = newFeH;
}

void Cluster::setMod (const double newMod)
{
    parameter[MOD] = newMod;
}

void Cluster::setAbs (const double newAbs)
{
    parameter[ABS] = newAbs;
}
