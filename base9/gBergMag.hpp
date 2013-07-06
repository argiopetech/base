#ifndef GBERGMAG_H
#define GBERGMAG_H

#include <string>

void loadBergeron (std::string path, MsFilter filterSet);
void bergeronTeffToMags (const std::vector<int> &filters, double wdLogTeff, double wdLogG, int wdType);

#endif
