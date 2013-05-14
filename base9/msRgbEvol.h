#if defined(MSRGBEVOL_H )
/* the file has already been defined */
#else
#define MSRGBEVOL_H

#include "evolve.h"
#include "msRgbEvol.h"
#include "gDsedMag.h"
#include "gGirMag.h"
#include "gChabMag.h"
#include "gYaleMag.h"
#include "structures.h"

void deriveAgbTipMass(struct cluster *pCluster);
double msRgbEvol(struct cluster *pCluster, double zamsMass);
double wdPrecLogAge(struct cluster *pCluster, double zamsMass);
void loadMSRgbModels(struct cluster *pCluster, char *path, int needFS);

#endif
