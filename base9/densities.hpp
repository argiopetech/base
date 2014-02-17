/* densities.h */

#ifndef DENSITIES_H
#define DENSITIES_H

#include <array>

#include "structures.hpp"
#include "Model.hpp"

const double DOF    = 6.0;
const double GAMMA6 = -2.0590305444197083417635;   /* GAMMA for DOF=6 */

double logPriorClust (const Cluster&, const Model&);
double logPost1Star (const Star&, const Cluster&, const Model&, const std::array<double, FILTS>&);
double logTDens (double, double, double, double);
#endif
