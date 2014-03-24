#include <array>
#include <iostream>
#include <fstream>
#include <mutex>
#include <string>
#include <sstream>
#include <thread>
#include <vector>
#include <stdexcept>

#include <boost/format.hpp>

#include "Star.hpp"
#include "densities.hpp"
#include "MsFilterSet.hpp"
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
vector<StellarSystem> readCmdData (struct ifmrMcmcControl &ctrl, const Model &evoModels, vector<int> &filters, std::array<double, FILTS> &filterPriorMin, std::array<double, FILTS> &filterPriorMax, const Settings &settings)
{
    string line, pch;

    vector<StellarSystem> systems;

    //Parse the header of the file to determine which filters are being used
    getline(ctrl.rData, line);  // Read in the header line

    istringstream header(line); // Ignore the first token (which is "id")

    header >> pch;

    while (!header.eof())
    {
        header >> pch;

        if (pch == "sig")
            break;                      // and check to see if they are 'sig'.  If they are, there are no more filters

        for (int filt = 0; filt < FILTS; filt++)
        {                               // Otherwise check to see what this filter's name is
            if (pch == evoModels.filterSet->getFilterName(filt))
            {
                filters.push_back(filt);
                const_cast<Model&>(evoModels).numFilts++;
                break;
            }
        }
    }

    if (filters.empty())
    {
        cerr << "Exiting due to empty filter set. Did you specify the correct filterSet?" << endl;
        exit(1);
    }

    // This loop reads in photometry data
    // It also reads a best guess for the mass
    systems.clear();

    while (getline(ctrl.rData, line))
    {
        systems.emplace_back(line, filters.size());

        for (decltype(filters.size()) i = 0; i < filters.size(); ++i)
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

    return systems;
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
