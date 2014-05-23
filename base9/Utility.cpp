#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "Utility.hpp"

using std::pair;
using std::ifstream;
using std::istringstream;
using std::string;
using std::vector;

namespace base
{
    namespace utility
    {
        void ThreadPool::parallelFor(const unsigned int size, std::function<void(const unsigned int)> func)
        {
            unsigned int loopThreads = (size < nThreads) ? size : nThreads;

            for (unsigned int idThread = 0; idThread < loopThreads; idThread++) {
                auto threadFunc = [=]() {
                    for (unsigned int i=idThread; i<size; i+=nThreads) {
                        func(i);
                    }
                };

                threads.at(idThread).addJob(threadFunc);
            }
            for (auto &t : threads)
                t.join();
        }

        // Assumptions:
        //   filterPriorMin and filterPriorMax will be cleared and overwritten.
        pair<vector<string>, vector<StellarSystem>> readPhotometry (ifstream &fin, vector<double> &filterPriorMin, vector<double> &filterPriorMax, const Settings &settings)
        {
            string line, pch;

            vector<StellarSystem> systems;
            vector<string> filterNames;

            //Parse the header of the file to determine which filters are being used
            getline(fin, line);  // Read in the header line

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
                throw InvalidPhotometry("Exiting due to empty filter set. Is your photometry file valid?");
            }

            // This loop reads in photometry data
            // It also reads a best guess for the mass
            systems.clear();

            /* Initialize filter prior mins and maxes */
            filterPriorMin.clear();
            filterPriorMin.resize(filterNames.size(), 1000);

            filterPriorMax.clear();
            filterPriorMax.resize(filterNames.size(), -1000);

            while (getline(fin, line))
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
        } // readPhotometry

    }
}
