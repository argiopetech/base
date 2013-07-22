#include <sstream>
#include <string>

#include <cmath>

#include "Cluster.hpp"
#include "Star.hpp"

using std::ifstream;
using std::string;
using std::stringstream;

double Star::getMass1() const
{
    return U;
}

double Star::getMass2() const
{
    return U * massRatio;
}

void Star::setMass1(double newMass)
{
    U = newMass;
}

// *** Unused ***
// void Star::setMass2 (const Cluster &pCluster, double newMass)
// {
//     massRatio = newMass / getMass1 (pCluster);
// }

void Star::readCMD(const string &s, int filters)
{
    double tempSigma;
    string starID;

    stringstream in(s);  

    in >> starID;

    for (int i = 0; i < filters; i++)
    {
        in >> obsPhot[i];
    }

    for (int i = 0; i < filters; i++)
    {
        in >> tempSigma;

        variance[i] = tempSigma * fabs (tempSigma);
        // The fabs() keeps the sign of the variance the same as that input by the user for sigma
        // Negative sigma (variance) is used to signal "don't count this band for this star"
    }

    in >> U
       >> massRatio
       >> status[0]
       >> clustStarPriorDens
       >> useDuringBurnIn;
}
