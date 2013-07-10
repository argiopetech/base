#ifndef GBERGMAG_H
#define GBERGMAG_H

#include <array>
#include <string>
#include <vector>

void loadBergeron (std::string path, MsFilter filterSet);
void bergeronTeffToMags (const std::vector<int> &filters, std::array<double, FILTS>&, double wdLogTeff, double wdLogG, int wdType);

#endif
