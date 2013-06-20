/* densities.h */

#ifdef DENSITIES_H
  /* the file has been included already */
#else
#define DENSITIES_H

#define DOF                  6.0
#define GAMMA6              -2.0590305444197083417635   /* GAMMA for DOF=6 */

double logPriorMass(struct star *p_Star, struct cluster *p_Clust);
double logPriorClust(struct cluster *p_Clust);
double logLikelihood(int numFilts, struct star *p_Star);
double tLogLikelihood(int numFilts, struct star *p_Star);
double scaledLogLike(int numFilts, struct star *pStar, double varScale);
double logPost1Star(struct star *p_Star, struct cluster *p_Clust);
double Phi(double x);
double logTDens(double x, double mean, double var, double nu);
#endif
