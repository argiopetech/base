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

void deriveAgbTipMass (struct cluster *pCluster);
double msRgbEvol (struct cluster *pCluster, double zamsMass);
double wdPrecLogAge (struct cluster *pCluster, double zamsMass);
void loadMSRgbModels (struct cluster *pCluster, std::string path, int needFS);

#endif
