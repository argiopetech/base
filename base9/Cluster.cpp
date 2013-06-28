
void readClust (FILE * pFile)
{
    int p;

    for (p = 0; p < NPARAMS; p++)
        fscanf (pFile, "%lg %lg %lg", &(pCluster->parameter[p]), &(pCluster->mean[p]), &(pCluster->stepSize[p]));

    fscanf (pFile, "%lf %lf %lf", &(pCluster->betamabs), &(pCluster->betaFabs), &(pCluster->betaFY));
}


void writeClust (FILE * pFile)
{
    int p;

    for (p = 0; p < NPARAMS; p++)
        fprintf (pFile, " %11.3e %11.3e %11.3e", pCluster->parameter[p], pCluster->mean[p], pCluster->stepSize[p]);
    fprintf (pFile, " %8.3f %8.3f %8.3f\n", pCluster->betamabs, pCluster->betaFabs, pCluster->betaFY);
}
