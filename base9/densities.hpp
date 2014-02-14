/* densities.h */

#ifndef DENSITIES_H
#define DENSITIES_H

#include <array>

#include "structures.hpp"
#include "Model.hpp"

const double DOF    = 6.0;
const double GAMMA6 = -2.0590305444197083417635;   /* GAMMA for DOF=6 */

double logPriorMass (const Star &p_Star, const Cluster &p_Clust);
double logPriorClust (const Cluster &p_Clust, const Model&);
double logLikelihood (int numFilts, const Star &p_Star, const std::array<double, FILTS>&, const std::array<double, FILTS>&);
double tLogLikelihood (int numFilts, const Star &p_Star, const std::array<double, FILTS>&, const std::array<double, FILTS>&);
double scaledLogLike (int numFilts, const Star &pStar, double varScale, const std::array<double, FILTS>&, const std::array<double, FILTS>&, const std::array<double, FILTS>&);
double logPost1Star (const Star &, const Cluster &, const Model&, const std::array<double, FILTS>&, const std::array<double, FILTS>&, const std::array<double, FILTS>&);
double Phi (double x);
double logTDens (double x, double mean, double var, double nu);
#endif
