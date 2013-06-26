#ifndef GBERGMAG_H
#define GBERGMAG_H

#include <string>

const int BERG_N_DA_LOG_G = 6;
const int BERG_N_DA_TEFF  = 57;
const int BERG_N_DB_LOG_G = 5;
const int BERG_N_DB_TEFF  = 31;
const int BERG_NFILTS     = 8;

void loadBergeron (std::string path, int filterSet);
void bergeronTeffToMags (double wdLogTeff, double wdLogG, int wdType);
#endif
