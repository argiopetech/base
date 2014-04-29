#include <array>
#include <iostream>
#include <string>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "Cluster.hpp"
#include "Star.hpp"

#include "DsedMsModel.hpp"
#include "LinearTransform.hpp"
#include "Matrix.hpp"

using std::array;
using std::string;
using std::vector;
using std::cerr;
using std::endl;


bool DsedMsModel::isSupported(FilterSetName filterSet) const
{
    return ( filterSet == FilterSetName::UBVRIJHK
          || filterSet == FilterSetName::SDSS
          || filterSet == FilterSetName::UVIS);
}


string DsedMsModel::getFileName (string path) const
{
    return path + "dsed/dsed_old.model";
}
