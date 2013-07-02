#ifndef DECIDE_H
#define DECIDE_H

#include <array>

#include "Model.hpp"

void decideFieldStar (Star stars1[], Cluster &pCluster, FILE * wFile, Model&);
void decideMass (Chain &mc, const Model&, std::array<double, 2>&);
void decideMassRatio (Chain &mc, const Model&, std::array<double, 2>&);
Cluster decideClust (Cluster clust1, Star stars1[], const int FS_ON_STATE, int *accept, int *reject, const int SAMPLE_TYPE, FILE * w_ptr, Model&, std::array<double, 2>&);
void updateVarScale (Star stars[], Cluster *pCluster);

#endif
