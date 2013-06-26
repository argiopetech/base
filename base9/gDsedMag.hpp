#ifndef GDSEDMAG_H
#define GDSEDMAG_H

#include <string>

const int N_DSED_Z         = 9;    /* number of metallicities in Chaboyer-Dotter isochrones */
const int N_DSED_AGES      = 52;   /* number of ages in Chaboyer-Dotter isochonres */
const int N_DSED_FILTS     = 8;
const int MAX_DSED_ENTRIES = 370;

void loadDsed (std::string path, int filterSet);
double deriveDsedAgbTip (double newFeH, double newY, double newAge);
double getDsedMags (double zamsMass);
double wdPrecLogAgeDsed (double thisFeH, double thisY, double zamsMass);
#endif
