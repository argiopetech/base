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

/*
 * Read data
 */
std::pair<vector<string>, vector<StellarSystem>> readCmdData (struct ifmrMcmcControl &ctrl, std::vector<double> &filterPriorMin, std::vector<double> &filterPriorMax, const Settings &settings)
{
    string line, pch;

    vector<StellarSystem> systems;
    vector<string> filterNames;

    //Parse the header of the file to determine which filters are being used
    getline(ctrl.rData, line);  // Read in the header line

    istringstream header(line); // Ignore the first token (which is "id")

    header >> pch;

    while (!header.eof())
    {
        header >> pch;

        // Check to see if input starts with "sig".  If yes, there are no more filters
        if (pch.compare(0, 3, "sig") == 0)
            break;

        filterNames.push_back(pch);
    }

    if (filterNames.empty())
    {
        cerr << "Exiting due to empty filter set. Did you specify the correct filterSet?" << endl;
        exit(1);
    }

    // This loop reads in photometry data
    // It also reads a best guess for the mass
    systems.clear();

    /* Initialize filter prior mins and maxes */
    filterPriorMin.resize(filterNames.size(), 1000);
    filterPriorMax.resize(filterNames.size(), -1000);

    while (getline(ctrl.rData, line))
    {
        systems.emplace_back(line, filterNames.size());

        for (size_t i = 0; i < filterNames.size(); ++i)
        {
            if (systems.back().obsPhot.at(i) < filterPriorMin.at(i))
            {
                filterPriorMin.at(i) = systems.back().obsPhot.at(i);
            }

            if (systems.back().obsPhot.at(i) > filterPriorMax.at(i))
            {
                filterPriorMax.at(i) = systems.back().obsPhot.at(i);
            }
        }

        if (!(systems.back().observedStatus == WD || (systems.back().obsPhot.at(settings.cluster.index) >= settings.cluster.minMag && systems.back().obsPhot.at(settings.cluster.index) <= settings.cluster.maxMag)))
        {
            systems.pop_back();
        }
    }

    return std::make_pair(filterNames, systems);
} /* readCmdData */


void printHeader (ofstream &file, array<double, NPARAMS> const &priors)
{
    const array<string, NPARAMS> paramNames = { "     logAge",
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
            file << paramNames.at(p) << ' ';
        }
    }
    file << "logPost" << endl;
}
