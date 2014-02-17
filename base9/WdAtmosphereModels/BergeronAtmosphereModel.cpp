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

void BergeronAtmosphereModel::loadModel (std::string path, MsFilter filterSet)
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
    double teff, ignore, mass;

    array<double, FILTS> mags;

    if (filterSet != MsFilter::UBVRIJHK && filterSet != MsFilter::SDSS)
    {
        cerr << "\nFilter set " << static_cast<int>(filterSet) << " not available on Bergeron models.  Exiting..." << endl;
        exit (1);
    }

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
            mags.fill(0.0);

            getline(fin, line);

            if (!line.empty())
            {
                if (line.at(2) != ' ') // This is an awful, kludgey way to check for the H/He split...
                {
                    stringstream in(line);

                    in >> teff
                       >> ignore >> ignore >> ignore;

                    for (int f = 0; f < 8; ++f) // Read the first 8 filters
                    {
                        in >> mags[f];
                    }

                    if (filterSet == MsFilter::SDSS) // iff we are using the ugrizJHK models
                    {
                        for (int f = 0; f < 5; ++f)
                        {
                            in >> mags[f]; // Read ugriz without overwriting JHK
                        }
                    }

                    // As long as we didn't run out of file somewhere in the middle, this doesn't trigger.
                    // Honestly, I think it should only happen if we have a corrupt file.
                    if (!fin.eof())
                    {
                        records.emplace_back(log10(teff), mags);
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

std::array<double, FILTS> BergeronAtmosphereModel::teffToMags (double wdLogTeff, double wdMass, WdAtmosphere wdType) const
{
    auto transformMags = [=](std::vector<struct AtmosCurve> curve)
        {
            Matrix<double, 2, BERG_NFILTS> logTeffMag;
            array<double, FILTS> mags;

            mags.fill(0.0);

            std::vector<AtmosCurve>::iterator massIter;

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

            std::array<std::vector<record>::iterator, 2> teffIter = {
                lower_bound(massIter[0].record.begin(), massIter[0].record.end(), record(wdLogTeff, mags)),
                lower_bound(massIter[1].record.begin(), massIter[1].record.end(), record(wdLogTeff, mags))
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
                for (int f = 0; f < BERG_NFILTS; ++f)
                {
                    logTeffMag[i][f] = linearTransform<>(teffIter[i][0].logTeff, teffIter[i][1].logTeff, teffIter[i][0].mags[f], teffIter[i][1].mags[f], wdLogTeff).val;
                }
            }

            //Interpolate in mass
            for (int f = 0; f < BERG_NFILTS; ++f)
            {
                mags[f] = linearTransform<>(massIter[0].mass, massIter[1].mass, logTeffMag[0][f], logTeffMag[1][f], wdMass).val;
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
