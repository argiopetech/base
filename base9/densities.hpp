/* densities.h */

#ifndef DENSITIES_H
#define DENSITIES_H

#include "structures.hpp"

const double DOF    = 6.0;
const double GAMMA6 = -2.0590305444197083417635;   /* GAMMA for DOF=6 */

double logPriorMass (Star *p_Star, Cluster *p_Clust);
double logPriorClust (Cluster *p_Clust);
double logLikelihood (int numFilts, Star *p_Star);
double tLogLikelihood (int numFilts, Star *p_Star);
double scaledLogLike (int numFilts, Star *pStar, double varScale);
double logPost1Star (Star *p_Star, Cluster *p_Clust);
double Phi (double x);
double logTDens (double x, double mean, double var, double nu);
#endif
