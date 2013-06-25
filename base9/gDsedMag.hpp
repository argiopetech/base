#if defined( GDSEDMAG_H )

#else
#define GDSEDMAG_H

//#define FILTS              8
#define N_DSED_Z           9    /* number of metallicities in Chaboyer-Dotter isochrones */
#define N_DSED_AGES        52   /* number of ages in Chaboyer-Dotter isochonres */
#define N_DSED_FILTS       8
#define MAX_DSED_ENTRIES   370

void loadDsed (char *path, int filterSet);
double deriveDsedAgbTip (double newFeH, double newY, double newAge);
double getDsedMags (double zamsMass);
double wdPrecLogAgeDsed (double thisFeH, double thisY, double zamsMass);
#endif
