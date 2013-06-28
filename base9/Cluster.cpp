#include <cstdio>

#include "Cluster.hpp"

void Cluster::readClust (FILE * pFile)
{
    int p;

    for (p = 0; p < NPARAMS; p++)
        fscanf (pFile, "%lg %lg %lg", &(parameter[p]), &(mean[p]), &(stepSize[p]));

    fscanf (pFile, "%lf %lf %lf", &(betamabs), &(betaFabs), &(betaFY));
}


void Cluster::writeClust (FILE * pFile)
{
    int p;

    for (p = 0; p < NPARAMS; p++)
        fprintf (pFile, " %11.3e %11.3e %11.3e", parameter[p], mean[p], stepSize[p]);

    fprintf (pFile, " %8.3f %8.3f %8.3f\n", betamabs, betaFabs, betaFY);
}
