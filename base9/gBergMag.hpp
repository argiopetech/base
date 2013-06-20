#if defined( GBERGMAG_H )

#else
#define GBERGMAG_H

//#define BERG_MASSES        11
//#define MAX_BERG_ENTRIES   100

void loadBergeron(char *path, int filterSet);
void bergeronTeffToMags(double wdLogTeff, double wdLogG, int wdType);

#define BERG_N_DA_LOG_G     6
#define BERG_N_DA_TEFF      57
#define BERG_N_DB_LOG_G     5
#define BERG_N_DB_TEFF      31
#define BERG_NFILTS         8

#endif
