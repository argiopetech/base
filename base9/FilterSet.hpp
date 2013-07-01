#include <array>

#include "constants.hpp"

void calcAbsCoeffs (MsFilterSet filterSet, std::array<double, FILTS> &clusterAbs);

void setFilterNames (MsFilterSet filterSet);
char *getFilterName (int index);
