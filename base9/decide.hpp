#ifndef DECIDE_H
#define DECIDE_H

void decideFieldStar (struct star stars1[], Cluster *pCluster, FILE * wFile);
void decideMass (struct chain *mc);
void decideMassRatio (struct chain *mc);
Cluster decideClust (Cluster clust1, struct star stars1[], const int FS_ON_STATE, int *accept, int *reject, const int SAMPLE_TYPE, FILE * w_ptr);
void updateVarScale (struct star stars[], Cluster *pCluster);

#endif
