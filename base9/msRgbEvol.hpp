#ifndef MSRGBEVOL_H
#define MSRGBEVOL_H

#include <string>

#include "evolve.hpp"
#include "msRgbEvol.hpp"
#include "gDsedMag.hpp"
#include "gGirMag.hpp"
#include "gChabMag.hpp"
#include "gYaleMag.hpp"
#include "structures.hpp"

void deriveAgbTipMass (Cluster *pCluster);
double msRgbEvol (Cluster *pCluster, double zamsMass);
double wdPrecLogAge (Cluster *pCluster, double zamsMass);
void loadMSRgbModels (Cluster *pCluster, std::string path, int needFS);

#endif
