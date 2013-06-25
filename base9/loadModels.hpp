#if defined( LOAD_MODELS_H )
/* the file has been included already */
#else
#define LOAD_MODELS_H

#include "msRgbEvol.h"
#include "gBergMag.h"
#include "wdCooling.h"
#include "gBaraffeMag.h"
#include "Settings.hpp"

void loadModels (int needFS, struct cluster *theCluster, struct Settings *);

#endif
