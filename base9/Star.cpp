#include "Cluster.hpp"
#include "Star.hpp"

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
