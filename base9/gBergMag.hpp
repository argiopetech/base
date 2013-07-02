#ifndef GBERGMAG_H
#define GBERGMAG_H

#include <string>

void loadBergeron (std::string path, MsFilterSet filterSet);
void bergeronTeffToMags (double wdLogTeff, double wdLogG, int wdType);

#endif
