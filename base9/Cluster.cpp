#include <cstdio>

#include "Cluster.hpp"
#include "structures.hpp"

double Cluster::getAge ()
{
    return parameter[AGE];
}

double Cluster::getY ()
{
    return parameter[YYY] + betaFY * (getFeH() - mean[FEH]);
}

double Cluster::getFeH ()
{
    return parameter[FEH];
}

double Cluster::getMod ()
{
    return parameter[MOD];
}

double Cluster::getAbs ()
{
    return parameter[ABS] + betaFabs * (getFeH() - mean[FEH]) + betamabs * (getMod() - mean[MOD]);
}

void Cluster::setAge (double newAge)
{
    parameter[AGE] = newAge;
}

void Cluster::setY (double newY)
{
    parameter[YYY] = newY - betaFY * (getFeH() - mean[FEH]);
}

void Cluster::setFeH (double newFeH)
{
    parameter[FEH] = newFeH;
}

void Cluster::setMod (double newMod)
{
    parameter[MOD] = newMod;
}

void Cluster::setAbs (double newAbs)
{
    parameter[ABS] = newAbs - betamabs * (getMod() - mean[MOD]) - betaFabs * (getFeH() - mean[FEH]);
}
