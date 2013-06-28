#ifndef DECIDE_H
#define DECIDE_H

void decideFieldStar (Star stars1[], Cluster *pCluster, FILE * wFile);
void decideMass (Chain *mc);
void decideMassRatio (Chain *mc);
Cluster decideClust (Cluster clust1, Star stars1[], const int FS_ON_STATE, int *accept, int *reject, const int SAMPLE_TYPE, FILE * w_ptr);
void updateVarScale (Star stars[], Cluster *pCluster);

#endif
