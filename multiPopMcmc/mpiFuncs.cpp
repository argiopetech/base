#include <array>
#include <iostream>
#include <fstream>
#include <mutex>
#include <string>
#include <sstream>
#include <thread>
#include <utility>
#include <vector>
#include <stdexcept>

#include <boost/format.hpp>

#include "Star.hpp"
#include "densities.hpp"
#include "Filters.hpp"
#include "marg.hpp"
#include "Model.hpp"
#include "mpiMcmc.hpp"
#include "samplers.hpp"
#include "Utility.hpp"
#include "WhiteDwarf.hpp"

using std::array;
using std::mutex;
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::isfinite;
using std::ofstream;
using std::istringstream;

using base::utility::ThreadPool;

void printHeader (ofstream &file, array<double, NPARAMS> const &priors)
{
    const array<string, NPARAMS> paramNames = {  "    logAge",
                                                "          Y",
                                                "        FeH",
                                                "    modulus",
                                                " absorption",
                                                "carbonicity",
                                                "  IFMRconst",
                                                "    IFMRlin",
                                                "   IFMRquad"};

    for (int p = 0; p < NPARAMS; p++)
    {
        if (priors.at(p) > EPSILON || p == MOD || p == FEH || p == ABS)
        {
            file << paramNames.at(p);
        }
    }
    file << "     logPost" << endl;
}
