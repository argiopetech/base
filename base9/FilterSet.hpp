#include <array>

#include "constants.hpp"

std::array<double, 8> calcAbsCoeffs (MsFilterSet filterSet);

void setFilterNames (MsFilterSet filterSet);
char *getFilterName (int index);
