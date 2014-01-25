#include <array>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>

#include <cmath>

#include "Cluster.hpp"
#include "Star.hpp"

#include "BergeronAtmosphereModel.hpp"
#include "evolve.hpp"
#include "LinearTransform.hpp"
#include "gBergMag.hpp"
#include "binSearch.hpp"

using std::array;
using std::ifstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

const unsigned int maxIgnore = std::numeric_limits<char>::max();

void BergeronAtmosphereModel::loadModel (std::string path, MsFilter filterSet)
{
    static struct std::string files[] = {
        "Table_Mass_0.2",
        "Table_Mass_0.3",
        "Table_Mass_0.4",
        "Table_Mass_0.5",
        "Table_Mass_0.6",
        "Table_Mass_0.7",
        "Table_Mass_0.8",
        "Table_Mass_0.9",
        "Table_Mass_1.0",
        "Table_Mass_1.2",
    };

    ifstream fin;
    string line, tempFile;
    double teff, logG;

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
        getline(fin, line);
        getline(fin, line);

        std::map<unsigned int, std::vector<struct record>> *theMap = &hCurves; // This lets us do the H and He curves in one loop

        while (!fin.eof())
        {
            // Teff logg Mbol BC U B V R I J H K u g r i z y b-y u-b v-y V-I G-R U-V U-G B-V Age
            double ignore;
            mags.fill(0.0);

            getline(fin, line);

            if (line.at(3) != ' ') // This is an awful, kludgey way to check for the H/He split...
            {
                stringstream in(line);

                in >> teff
                   >> logG
                   >> ignore >> ignore;

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
                    int tTeff = static_cast<int>(teff);
                    (*theMap)[tTeff].emplace_back(logG, mags);
                }
            }
            else // This is the split point between H and He tables
            {
                getline(fin, line); // Eat the extra header line
                theMap = &heCurves; // And swap maps to the helium curve
            }
        }

        fin.close();
    }
}

std::array<double, FILTS> BergeronAtmosphereModel::teffToMags (double wdLogTeff, double wdLogG, WdAtmosphere wdType)
{
    auto transformMags = [](std::map<unsigned int, std::vector<struct record>> theMap)
        {
            Matrix<2, BERG_NFILTS> logGMag;

            array<double, 2> logGMag;
            array<double, FILTS> mags;

            mags.fill(0.0);

            if (wdLogG < bLogG[0])
                l = 0;
            else if (wdLogG >= bLogG[BERG_N_DA_LOG_G - 2])
                l = BERG_N_DA_LOG_G - 2;
            else
                l = (int) floor ((wdLogG - bLogG[0]) / (bLogG[1] - bLogG[0]));

            auto t = theMap.lower_bound(wdLogTeff);

            if (teffIter == theMap.begin()) {
                // log << "Poetential Teff underflow in WD Atmosphere Model" << endl;
            }
            else if (teffIter == theMap.end()) {
                // log << "Teff overflow in WD Atmosphere Model" << endl;
                teffIter -= 2;
            }
            else {
                teffIter -= 1;
            }

            auto gIter = { lower_bound(teffIter[0].second.begin(), teffIter[0].second.end(), wdLogG),
                           lower_bound(teffIter[1].second.begin(), teffIter[1].second.end(), wdLogG)};

            for ( int i = 0; i < 2; ++i )
            {
                if (gIter[i] == teffIter[i].second.begin()) {
                    // log << "Poetential Teff underflow in WD Atmosphere Model" << endl;
                }
                else if (gIter[i] == teffIter[i].second.end()) {
                    // log << "Teff overflow in WD Atmosphere Model" << endl;
                    gIter[i] -= 2;
                }
                else {
                    gIter[i] -= 1;
                }
            }

            //Interpolate in logTeff
            for (int i = 0; i < 2; i++)
            {
                logGTrans[i] = linearTransform<>(teffIter[0].first, teffIter[1].first, gIter[i].logG, gIter[i + 1].logG, wdLogTeff).val;

                for (int f = 0; f < BERG_NFILTS; ++f)
                {
                    logGMag[i][f] = linearTransform<>(teffIter[0].first, teffIter[1].first, gIter[i].mags[f], gIter[i + 1].mags[f], wdLogTeff).val;
                }
            }

            //Interpolate in log(g)
            for (int f = 0; f < BERG_NFILTS; ++f)
            {
                mags[f] = linearTransform<>(logGTrans[0], logGTrans[1], logGMag[0][f], logGMag[1][f], wdLogG).val;
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
        std::cerr << "Somehow, you've input an invalide WdAtmosphere type." << std::endl;
        exit(-1); // I know it isn't nice to exit, but this should be possible in the first place.
    }
}
