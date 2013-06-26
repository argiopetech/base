/* densities.h */

#ifndef DENSITIES_H
#define DENSITIES_H

const double DOF    = 6.0;
const double GAMMA6 = -2.0590305444197083417635;   /* GAMMA for DOF=6 */

double logPriorMass (struct star *p_Star, struct cluster *p_Clust);
double logPriorClust (struct cluster *p_Clust);
double logLikelihood (int numFilts, struct star *p_Star);
double tLogLikelihood (int numFilts, struct star *p_Star);
double scaledLogLike (int numFilts, struct star *pStar, double varScale);
double logPost1Star (struct star *p_Star, struct cluster *p_Clust);
double Phi (double x);
double logTDens (double x, double mean, double var, double nu);
#endif
