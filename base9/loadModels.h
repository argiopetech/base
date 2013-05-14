#if defined( LOAD_MODELS_H )
  /* the file has been included already */
#else
#define LOAD_MODELS_H

#include "msRgbEvol.h"
#include "gBergMag.h"
#include "wdCooling.h"
#include "gBaraffeMag.h"

void loadModels(int needFS, struct cluster *theCluster);

#endif
