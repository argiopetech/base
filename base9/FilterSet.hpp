#include <array>

#include "constants.hpp"

std::array<double, FILTS> calcAbsCoeffs (MsFilterSet filterSet);

void setFilterNames (MsFilterSet filterSet);
char *getFilterName (int index);
