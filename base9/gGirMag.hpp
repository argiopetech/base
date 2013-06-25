#if defined( GGIRMAG_H )

#else
#define GGIRMAG_H

#define N_GIR_Z            8    /* number of metallicities in Girardi isochrones */
#define N_GIR_AGES         50   /* number of ages in Girardi isochonres */
#define N_GIR_FILTS        8
#define MAX_GIR_ENTRIES    10000        /* max Girardi entries for given Z */

void loadGirardi (char *path, int filterSet);
double getGirardiMags (double zamsMass);
double interpInMass (int a, double zamsMass, int indexFeH, double *tempMag);
double deriveGirAgbTip (double newFeH, double newY, double newAge);
double wdPrecLogAgeGir (double thisFeH, double thisY, double zamsMass);

#endif
