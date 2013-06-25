#if defined( GCHABMAG_H )

#else
#define GCHABMAG_H

#define N_CHAB_FILTS       8
#define N_CHAB_Z           4    /* number of metallicities in Chaboyer-Dotter isochrones */
#define N_CHAB_Y           5    /* number of He abundances in Chaboyer-Dotter isochrones */
#define N_CHAB_AGES        19   /* number of ages in Chaboyer-Dotter isochonres */
#define LOW                0
#define HIGH               1
#define MAX_CHAB_ENTRIES   280

void loadChaboyer (char *path, int filterSet);
double deriveChabAgbTip (double newFeH, double newY, double newAge);
double getChaboyerMags (double zamsMass);
double wdPrecLogAgeChaboyer (double FeH, double thisY, double zamsMass);
#endif
