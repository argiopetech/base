#include <sstream>
#include <string>

#include <cmath>

#include "Cluster.hpp"
#include "Star.hpp"

using std::ifstream;
using std::string;
using std::stringstream;

double Star::getMass1 (const Cluster &pCluster) const
{
    return U + beta[AGE][0] * (pCluster.getAge() - pCluster.mean[AGE])
             + beta[MOD][0] * (pCluster.getMod() - pCluster.mean[MOD])
             + beta[FEH][0] * (pCluster.getFeH() - pCluster.mean[FEH])
             + beta[YYY][0] * (pCluster.getY()   - pCluster.mean[YYY])
             + betaMassRatio[0] * pow (massRatio, betaMassRatio[1]
             );
}

double Star::getMass2 (const Cluster &pCluster) const
{
    return getMass1 (pCluster) * massRatio;
}

void Star::setMass1 (const Cluster &pCluster, double newMass)
{
    U = newMass - ( beta[AGE][0] * (pCluster.getAge() - pCluster.mean[AGE]) 
                  + beta[MOD][0] * (pCluster.getMod() - pCluster.mean[MOD]) 
                  + beta[FEH][0] * (pCluster.getFeH() - pCluster.mean[FEH]) 
                  + beta[YYY][0] * (pCluster.getY()   - pCluster.mean[YYY])
                  + betaMassRatio[0] * pow (massRatio, betaMassRatio[1])
                  );
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
