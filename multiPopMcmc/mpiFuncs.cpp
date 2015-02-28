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

#include "Star.hpp"
#include "densities.hpp"
#include "Filters.hpp"
#include "marg.hpp"
#include "Model.hpp"
#include "mpiMcmc.hpp"
#include "samplers.hpp"
#include "Utility.hpp"
#include "WhiteDwarf.hpp"
#include "MpiMcmcApplication.hpp"

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
    const array<string, NPARAMS> paramNames = {  "     logAge",
                                                "          YA",
                                                "         FeH",
                                                "     modulus",
                                                "  absorption",
                                                " carbonicity",
                                                "          YB",
                                                "      Lambda",
                                                "    IFMRquad"};

    for (int p = 0; p < NPARAMS; p++)
    {
        if ((priors.at(p) > EPSILON || p == MOD || p == FEH || p == ABS) && !(p == YYA || p == YYB || p == LAMBDA))
        {
            file << paramNames.at(p);
        }
    }
    file << paramNames.at(YYA) << paramNames.at(YYB) << paramNames.at(LAMBDA);
    file << "     logPost" << endl;
}
