#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "sqlite/sqlite3.h"
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
        namespace sqlite
        {
            void prepQuery(sqlite3 *db, sqlite3_stmt **stmt, string sql)
            {
                while ( SQLITE_BUSY == sqlite3_prepare_v2( db, sql.c_str(), -1,
                                                           stmt, nullptr ))
                { ; }
            }

            int stepQuery(sqlite3_stmt *stmt)
            {
                int ret;

                do
                {
                    ret = sqlite3_step(stmt);
                } while ( ret == SQLITE_BUSY );

                return ret;
            }

            void finalize(sqlite3_stmt *stmt)
            {
                sqlite3_finalize(stmt);
                stmt = nullptr;
            }

            void dbError(sqlite3 *db)
            {
                throw std::runtime_error("SQLite failed: " + string(sqlite3_errmsg(db)));
            }
        }

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


        // Assumptions:
        //   filterPriorMin and filterPriorMax will be cleared and overwritten.
        pair<vector<string>, vector<StellarSystem>> readPhotometryFromDB (vector<double> &filterPriorMin, vector<double> &filterPriorMax, const Settings &settings)
        {
            using namespace sqlite;

            string line, pch;

            vector<StellarSystem> systems;
            vector<string> filterNames;

            sqlite3 *db;
            sqlite3_stmt* stmt;
            sqlite3_stmt* photometry;

            string dbName = settings.files.output + ".base9";

            if (sqlite3_open(dbName.c_str(), &db) != SQLITE_OK)
            {
                throw std::runtime_error("Failed to open database " + dbName + ".\nSQLite says: " + sqlite3_errmsg(db));
            }

            //Run a query to determine which filters are being used
            {
                prepQuery(db, &stmt,
                          "select filter from photometry where runid = ?1 group by filter order by filter;");
                sqlite3_bind_int(stmt, 1, settings.run);

                int ret = stepQuery(stmt);

                while (ret == SQLITE_ROW)
                {
                    filterNames.emplace_back(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0)));

                    ret = stepQuery(stmt);
                }

                if (ret != SQLITE_DONE)
                {
                    dbError(db);
                }

                finalize(stmt);
            }

            if (filterNames.empty())
            {
                throw InvalidPhotometry("Exiting due to empty filter set. Is this a valid run?");
            }

            // This loop reads in photometry data
            // It also reads a best guess for the mass
            systems.clear();

            /* Initialize filter prior mins and maxes */
            filterPriorMin.clear();
            filterPriorMin.resize(filterNames.size(), 1000);

            filterPriorMax.clear();
            filterPriorMax.resize(filterNames.size(), -1000);

            {
                // Prepare query to get stars
                prepQuery(db, &stmt,
                          "select starid, primaryMass, secondaryMassRatio, stage, clusterMembershipPrior, "
                          "useDuringBurnIn from star where runid = ?1;");
                sqlite3_bind_int(stmt, 1, settings.run);

                // Prepare query to get filters
                prepQuery(db, &photometry,
                          "select magnitude, stdDev from photometry "
                          "where runid = ?1 and starid = ?2 order by filter;");

                int ret = stepQuery(stmt);

                while (ret == SQLITE_ROW)
                {
                    vector<double> obsPhot;
                    vector<double> stdDevs;

                    string id = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
                    double mass       = sqlite3_column_double(stmt, 1);
                    double massRatio  = sqlite3_column_double(stmt, 2);
                    int    status     = sqlite3_column_int   (stmt, 3);
                    double clustPrior = sqlite3_column_double(stmt, 4);
                    bool   useDBI     = sqlite3_column_int   (stmt, 5);

                    // Run query for obsPhot and stdDevs;
                    sqlite3_reset(photometry);
                    sqlite3_bind_int(photometry, 1, settings.run);
                    sqlite3_bind_text(photometry, 2, id.data(), id.length(), SQLITE_STATIC);

                    int photRet = stepQuery(photometry);

                    while (photRet == SQLITE_ROW)
                    {
                        obsPhot.push_back(sqlite3_column_double(photometry, 0));
                        stdDevs.push_back(sqlite3_column_double(photometry, 1));

                        photRet = stepQuery(photometry);
                    }

                    if (photRet != SQLITE_DONE)
                    {
                        dbError(db);
                    }

                    StellarSystem system;
                    system.setSystemParams(id, obsPhot, stdDevs, mass, massRatio,
                                           status, clustPrior, useDBI);

                    systems.push_back(system);

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
                             ||(systems.back().obsPhot.at(settings.cluster.index) >= settings.cluster.minMag
                               && systems.back().obsPhot.at(settings.cluster.index) <= settings.cluster.maxMag)))
                        {
                            systems.pop_back();
                        }
                    }
                    catch (std::out_of_range &e)
                    {
                        throw std::out_of_range(string(e.what()) + " while loading photometry. Check your filter index.");
                    }

                    ret = stepQuery(stmt);
                }

                if (ret != SQLITE_DONE)
                {
                    dbError(db);
                }
            }

            sqlite3_close(db);
            db = nullptr;

            return std::make_pair(filterNames, systems);
        } // readPhotometryFromDB

        vector<SinglePopMcmcRecord> readSinglePopMainRunResFromDB(int run, const string filepath)
        {
            using namespace sqlite;

            vector<SinglePopMcmcRecord> results;

            sqlite3 *db;
            sqlite3_stmt* stmt;

            string dbName = filepath + ".base9";

            if (sqlite3_open(dbName.c_str(), &db) != SQLITE_OK)
            {
                throw std::runtime_error("Failed to open database " + dbName + ".\nSQLite says: " + sqlite3_errmsg(db));
            }

            prepQuery(db, &stmt,
                      "select iterId, logAge, y, feh, modulus, absorption, carbonFraction, "
                      "ifmrIntercept, ifmrSlope, ifmrQuadCoef, logPost, stage "
                      "from single_pop "
                      "where runid = ?1 and stage >= 2;");
            sqlite3_bind_int(stmt, 1, run);

            int ret = stepQuery(stmt);

            while (ret == SQLITE_ROW)
            {
                Cluster clust;

                auto iter = sqlite3_column_int   (stmt, 0);
                clust.age = sqlite3_column_double(stmt, 1);
                clust.yyy = sqlite3_column_double(stmt, 2);
                clust.feh = sqlite3_column_double(stmt, 3);
                clust.mod = sqlite3_column_double(stmt, 4);
                clust.abs = sqlite3_column_double(stmt, 5);

                clust.carbonicity   = sqlite3_column_double(stmt, 6);
                clust.ifmrIntercept = sqlite3_column_double(stmt, 7);
                clust.ifmrSlope     = sqlite3_column_double(stmt, 8);
                clust.ifmrQuadCoef  = sqlite3_column_double(stmt, 9);

                auto logP = sqlite3_column_double(stmt, 10);
                auto stag = sqlite3_column_double(stmt, 11);

                results.push_back({{iter}, static_cast<AdaptiveMcmcStage>(stag), clust, logP, false});

                ret = stepQuery(stmt);
            }

            if (ret != SQLITE_DONE)
            {
                dbError(db);
            }

            finalize(stmt);

            sqlite3_close(db);
            db = nullptr;

            return results;
        }

        // This sets up formatting suitable for either a string or a double in the default column output format
        // This is the least efficient possible way of doing this, but output shouldn't be a bottlneck
        std::ostream& format (std::ostream& out)
        {
            return out << std::setw(11) << std::fixed << std::setprecision(6);
        }
    }
}
