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
using std::stringstream;
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
                    for (unsigned int i = idThread; i < size; i += nThreads) {
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

                try
                {
                    if (!(systems.back().observedStatus == WD
                       ||   (systems.back().obsPhot.at(settings.cluster.index) >= settings.cluster.minMag
                          && systems.back().obsPhot.at(settings.cluster.index) <= settings.cluster.maxMag)))
                    {
                        systems.pop_back();
                    }
                }
                catch (std::out_of_range &e)
                {
                    throw std::out_of_range(string(e.what()) + " while loading photometry. Check your filter index.");
                }
            }

            return std::make_pair(filterNames, systems);
        } // readPhotometry

        // This sets up formatting suitable for either a string or a double in the default column output format
        // This is the least efficient possible way of doing this, but output shouldn't be a bottlneck
        std::ostream& format (std::ostream& out)
        {
            return out << std::setw(11) << std::fixed << std::setprecision(6);
        }


        /*
         * Read sampled params
         */
        vector<clustPar> readSampledParams (Model &evoModels, const Settings &s)
        {
            auto acceptStage = [&s](int stage){ return (stage == 3) || s.includeBurnin; };

            return readSampledParams(s.files.output + ".res", evoModels, s, acceptStage);
        }

        vector<clustPar> readSampledParams (string filename, Model &evoModels, const Settings &s, std::function<bool(int)> acceptStage)
        {
            vector<clustPar> sampledPars;

            string line;

            std::ifstream parsFile;
            parsFile.open(filename);

            bool hasY, hasCarbonicity;

            getline(parsFile, line); // Parse header

            {
                string sin;
                stringstream in(line);

                in >> sin  // logAge
                   >> sin; // Y?

                if (sin == "Y")
                {
                    hasY = true;

                    in >> sin  // FeH
                       >> sin  // Abs
                       >> sin  // DistMod
                       >> sin; // Carbonicity?

                    if (sin == "carbonicity")
                        hasCarbonicity = true;
                    else
                        hasCarbonicity = false;
                }
                else
                {
                    hasY = false;

                    // This one skips FeH (because it was already read instead of Y)
                    in >> sin  // Abs
                       >> sin  // DistMod
                       >> sin; // Carbonicity?

                    if (sin == "carbonicity")
                        hasCarbonicity = true;
                    else
                        hasCarbonicity = false;
                }
            }

            while (getline(parsFile, line))
            {
                stringstream in(line);

                double newAge, newY, newFeh, newMod, newAbs, newCarbonicity, newIInter, newISlope, newIQuad, newLogPost;
                int stage;

                in >> newAge;

                if (hasY)
                    in >> newY;
                else
                    newY = s.cluster.starting.Y;

                in >> newFeh
                   >> newMod
                   >> newAbs;

                if (hasCarbonicity)
                    in >> newCarbonicity;
                else
                    newCarbonicity = s.cluster.starting.carbonicity;

                if (evoModels.IFMR >= 4 && evoModels.IFMR < 12)
                {
                    in >> newIInter
                       >> newISlope;
                }

                if (evoModels.IFMR >= 9 && evoModels.IFMR < 12)
                {
                    in >> newIQuad;
                }

                in >> newLogPost;

                in >> stage;

                auto amsStage = (AdaptiveMcmcStage)stage;

                if (acceptStage(stage))
                {
                    sampledPars.emplace_back(newAge, newY, newFeh, newMod, newAbs, newCarbonicity, newIInter, newISlope, newIQuad, newLogPost, amsStage);
                }
            }

            parsFile.close();

            return sampledPars;
        }
    }
}
