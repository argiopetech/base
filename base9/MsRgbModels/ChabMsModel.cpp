#include <string>
#include <iostream>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "Cluster.hpp"
#include "Star.hpp"

#include "ChabMsModel.hpp"
#include "LinearTransform.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;

bool ChabMsModel::isSupported(FilterSetName filterSet) const
{
    return (filterSet == FilterSetName::UBVRIJHK);
}

string ChabMsModel::getFileName (string path) const
{
    return path + "chaboyer/chaboyer.model";
}
