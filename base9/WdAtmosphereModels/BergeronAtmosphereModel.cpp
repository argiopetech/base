#include <array>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include <cmath>

#include "Cluster.hpp"
#include "Star.hpp"

#include "BergeronAtmosphereModel.hpp"
#include "LinearTransform.hpp"
#include "Matrix.hpp"

using std::array;
using std::ifstream;
using std::stringstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

void BergeronAtmosphereModel::loadModel (std::string path)
{
    static std::string files[] = {
        {"Table_Mass_0.2"},
        {"Table_Mass_0.3"},
        {"Table_Mass_0.4"},
        {"Table_Mass_0.5"},
        {"Table_Mass_0.6"},
        {"Table_Mass_0.7"},
        {"Table_Mass_0.8"},
        {"Table_Mass_0.9"},
        {"Table_Mass_1.0"},
        {"Table_Mass_1.2"}
    };

    ifstream fin;
    string line, tempFile;
    double teff, ignore, mass, logg;

    // Open the appropriate file for each mass
    for (auto f : files)
    {
        tempFile = path + "bergeron/" + f;

        fin.open(tempFile);

        if (!fin)
        {
            cerr << "\nFile " << tempFile << " was not found - exiting" << endl;
            exit (1);
        }

        //Skip header lines
        fin >> ignore >> mass; // Get the mass out of line 1 (column 2)

        getline(fin, line); // Eat the rest of the first line
        getline(fin, line);

        std::vector<record> records;

        while (!fin.eof())
        {
            // Teff logg Mbol BC U B V R I J H K u g r i z y b-y u-b v-y V-I G-R U-V U-G B-V Age
            vector<double> mags;

            getline(fin, line);

            if (!line.empty())
            {
                if (line.at(2) != ' ') // This is an awful, kludgey way to check for the H/He split...
                {
                    stringstream in(line);

                    in >> teff
                       >> logg
                       >> ignore >> ignore;

                    for (auto f : availableFilters) // Read all filters
                    {
                        double tmp;
                        in >> tmp;
                        mags.push_back(tmp);
                    }

                    // As long as we didn't run out of file somewhere in the middle, this doesn't trigger.
                    // Honestly, I think it should only happen if we have a corrupt file.
                    if (!fin.eof())
                    {
                        records.emplace_back(log10(teff), logg, mags);
                    }
                }
                else // This is the split point between H and He tables
                {
                    hCurves.emplace_back(mass, records);
                    records.clear();

                    getline(fin, line); // Eat the extra header line
                }
            }
        }

        heCurves.emplace_back(mass, records);

        fin.close();
    }
}

std::vector<double> BergeronAtmosphereModel::teffToMags (const double wdLogTeff, double wdMass, WdAtmosphere wdType) const
{
    auto transformMags = [=](const std::vector<struct AtmosCurve> &curve)
        {
            MVatrix<double, 2> logTeffMag;
            vector<double> mags;

            vector<AtmosCurve>::const_iterator massIter;

            if (wdMass < 0.2) // Below the hard-coded bottom mass
            {
                massIter = curve.begin();
            }
            else if (wdMass >= 1.0) // Above the hard-coded top - 2 (there is no 1.1 Mo curve)
            {
                massIter = curve.end() - 2;
            }
            else // Otherwise, hard-coded math magic.
            {
                massIter = curve.begin() + (static_cast<int>(wdMass * 10) - 2);
            }

            array<vector<record>::const_iterator, 2> teffIter = {
                lower_bound(massIter[0].record.begin(), massIter[0].record.end(), record(wdLogTeff, 0, mags)),
                lower_bound(massIter[1].record.begin(), massIter[1].record.end(), record(wdLogTeff, 0, mags))
            };

            for ( int i = 0; i < 2; ++i )
            {
                if (teffIter[i] == massIter[i].record.begin()) {
                    // log << "Poetential Teff underflow in WD Atmosphere Model" << endl;
                }
                else if (teffIter[i] == massIter[i].record.end()) {
                    // log << "Teff overflow in WD Atmosphere Model" << endl;
                    teffIter[i] -= 2;
                }
                else {
                    teffIter[i] -= 1;
                }
            }

            //Interpolate in logTeff
            for (int i = 0; i < 2; i++)
            {
                for (size_t f = 0; f < teffIter[i][0].mags.size(); ++f)
                {
                    logTeffMag[i].push_back(linearTransform<>(teffIter[i][0].logTeff, teffIter[i][1].logTeff, teffIter[i][0].mags[f], teffIter[i][1].mags[f], wdLogTeff).val);
                }
            }

            //Interpolate in mass
            for (size_t f = 0; f < logTeffMag[0].size(); ++f)
            {
                mags.push_back(linearTransform<>(massIter[0].mass, massIter[1].mass, logTeffMag[0][f], logTeffMag[1][f], wdMass).val);
            }

            return mags;
        };

    if (wdType == WdAtmosphere::DA)
    {
        return transformMags(hCurves);
    }
    else if (wdType == WdAtmosphere::DB)
    {
        return transformMags(heCurves);
    }
    else
    {
        std::cerr << "Somehow, you've input an invalid WdAtmosphere type." << std::endl;
        exit(-1); // I know it isn't nice to exit, but this should be possible in the first place.
    }
}

double BergeronAtmosphereModel::teffToLogg (double wdLogTeff, double wdMass, WdAtmosphere wdType) const
{
    auto transformMags = [=](const std::vector<struct AtmosCurve> &curve)
        {
            array<double, 2> logTeffMag;

            vector<AtmosCurve>::const_iterator massIter;

            if (wdMass < 0.2) // Below the hard-coded bottom mass
            {
                massIter = curve.begin();
            }
            else if (wdMass >= 1.0) // Above the hard-coded top - 2 (there is no 1.1 Mo curve)
            {
                massIter = curve.end() - 2;
            }
            else // Otherwise, hard-coded math magic.
            {
                massIter = curve.cbegin() + (static_cast<int>(wdMass * 10) - 2);
            }

            array<vector<record>::const_iterator, 2> teffIter = {
                lower_bound(massIter[0].record.begin(), massIter[0].record.end(), record(wdLogTeff, 0, {})),
                lower_bound(massIter[1].record.begin(), massIter[1].record.end(), record(wdLogTeff, 0, {})),
            };

            for ( int i = 0; i < 2; ++i )
            {
                if (teffIter[i] == massIter[i].record.begin()) {
                    // log << "Poetential Teff underflow in WD Atmosphere Model" << endl;
                }
                else if (teffIter[i] == massIter[i].record.end()) {
                    // log << "Teff overflow in WD Atmosphere Model" << endl;
                    teffIter[i] -= 2;
                }
                else {
                    teffIter[i] -= 1;
                }
            }

            //Interpolate in logTeff
            for (int i = 0; i < 2; i++)
            {
                logTeffMag[i] = linearTransform<>(teffIter[i][0].logTeff
                                                , teffIter[i][1].logTeff
                                                , teffIter[i][0].logg
                                                , teffIter[i][1].logg
                                                , wdLogTeff).val;
            }

            return linearTransform<>(massIter[0].mass
                                   , massIter[1].mass
                                   , logTeffMag[0]
                                   , logTeffMag[1]
                                   , wdMass).val;
        };

    if (wdType == WdAtmosphere::DA)
    {
        return transformMags(hCurves);
    }
    else if (wdType == WdAtmosphere::DB)
    {
        return transformMags(heCurves);
    }
    else
    {
        std::cerr << "Somehow, you've input an invalid WdAtmosphere type." << std::endl;
        exit(-1); // I know it isn't nice to exit, but this should be possible in the first place.
    }
}

void BergeronAtmosphereModel::restrictToFilters(const vector<string>& filters)
{
    vector<int> indices;

    for (auto f : filters)
    {
        bool foundFilter = false;

        for (size_t i = 0; i < availableFilters.size(); ++i)
        {
            if ( f == availableFilters.at(i) )
            {
                indices.push_back(i);

                foundFilter = true;
                break;
            }
        }

        if ( ! foundFilter )
        {
            throw InvalidModelError(f);
        }
    }

    for (size_t i = 0; i < indices.size(); ++i)
    {
        // Prepend the desired filter values to the beginning of the list...
        availableFilters.insert(availableFilters.begin() + i, availableFilters.at(indices.at(i) + i));
    }

    // ...then erase the rest.
    availableFilters.erase(availableFilters.begin() + indices.size(), availableFilters.end());
    // And try to reduce the capacity() to what is actually needed
    availableFilters.shrink_to_fit();

    // We should now have the same number of filters in the available filter set as in the input vector.
    assert(availableFilters.size() == filters.size());

    // Having this as a seperate loop allows a size assertion without
    // going over the data structure again. It also probably isn't bad
    // for cache locality (not that it matters)
    for (auto &atmos : hCurves)
    {
        for (auto &r : atmos.record)
        {
            for (size_t i = 0; i < indices.size(); ++i)
            {
                // Prepend the desired filter values to the beginning of the list...
                r.mags.insert(r.mags.begin() + i, r.mags.at(indices.at(i) + i));
            }

            // ...then erase the rest.
            r.mags.erase(r.mags.begin() + indices.size(), r.mags.end());
            // And try to reduce the capacity() to what is actually needed
            r.mags.shrink_to_fit();

            // At this point, eep.mags, availableFilters, and filters should all be the same size.
            assert(r.mags.size() <= filters.size());
        }
    }

    for (auto &atmos : heCurves)
    {
        for (auto &r : atmos.record)
        {
            for (size_t i = 0; i < indices.size(); ++i)
            {
                // Prepend the desired filter values to the beginning of the list...
                r.mags.insert(r.mags.begin() + i, r.mags.at(indices.at(i) + i));
            }

            // ...then erase the rest.
            r.mags.erase(r.mags.begin() + indices.size(), r.mags.end());
            // And try to reduce the capacity() to what is actually needed
            r.mags.shrink_to_fit();

            // At this point, eep.mags, availableFilters, and filters should all be the same size.
            assert(r.mags.size() <= filters.size());
        }
    }
}
