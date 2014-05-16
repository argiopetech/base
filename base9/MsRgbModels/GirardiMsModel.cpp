#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "Cluster.hpp"
#include "Isochrone.hpp"
#include "Star.hpp"

#include "LinearTransform.hpp"
#include "GirardiMsModel.hpp"

using std::array;
using std::string;
using std::vector;
using std::stringstream;
using std::lower_bound;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;

string GirardiMsModel::getFileName (string path) const
{
    return path + "girardi/girardi.model";
}


bool GirardiMsModel::isSupported(FilterSetName filterSet) const
{
    return (filterSet == FilterSetName::UBVRIJHK || filterSet == FilterSetName::ACS);
}
