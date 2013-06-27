#ifndef DECIDE_H
#define DECIDE_H

void decideFieldStar (struct star stars1[], struct cluster *pCluster, FILE * wFile);
void decideMass (struct chain *mc);
void decideMassRatio (struct chain *mc);
struct cluster decideClust (struct cluster clust1, struct star stars1[], const int FS_ON_STATE, int *accept, int *reject, const int SAMPLE_TYPE, FILE * w_ptr);
void updateVarScale (struct star stars[], struct cluster *pCluster);

#endif
