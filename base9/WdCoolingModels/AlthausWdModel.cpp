#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>

#include "Cluster.hpp"
#include "Star.hpp"

#include "AlthausWdModel.hpp"

using std::lower_bound;
using std::string;
using std::ifstream;
using std::vector;
using std::cerr;
using std::endl;

struct althausModel
{
    const string filename;
    const bool hasHLum;
    const double mass;
};

void AlthausWdModel::loadModel(string path, FilterSetName)
{
    static struct althausModel althaus[] = {
        {"T045_1E4.Z0", false, 0.45},
        {"T047_1E4.Z0", false, 0.47},
        {"T05_1E4.Z0", false, 0.50},
        {"T052_1E4.Z0", false, 0.52},
        {"T054_1E4.Z0", false, 0.54},
        {"T056_1E4.Z0", false, 0.56},
        {"T058_1E4.Z0", true, 0.58},
        {"T06_1E4.Z0", true, 0.60},
        {"T062_1E4.Z0", true, 0.62},
        {"T064_1E4.Z0", true, 0.64},
        {"T066_1E4.Z0", true, 0.66},
        {"T068_1E4.Z0", true, 0.68},
        {"T07_1E4.Z0", true, 0.70},
        {"T072_1E4.Z0", true, 0.72},
        {"T074_1E4.Z0", true, 0.74},
        {"T076_1E4.Z0", true, 0.76},
        {"T078_1E4.Z0", true, 0.78},
        {"T08_1E4.Z0", true, 0.80},
        {"T082_1E4.Z0", true, 0.82},
        {"T084_1E4.Z0", true, 0.84},
        {"T09_1E4.Z0", true, 0.90},
        {"T10_1E4.Z0", true, 1.00},
        {"T11_1E4.Z0", true, 1.10},
    };

    string tempFile, line;
    double newAge, newTeff, newRadius;
    double newCarbon = 0.6 ;// 0.38; // Good default value, per Mike Montgomery
    double ignore;

    for (auto massModel : althaus)
    {
        tempFile = path + "althaus/" + massModel.filename;

        ifstream pCoolingModels;
        pCoolingModels.open(tempFile);

        if (!pCoolingModels.is_open())
        {
            cerr << "\n file " << tempFile << " was not found - exiting" << endl;
            exit (1);
        }        

        getline(pCoolingModels, line); // get header line
        getline(pCoolingModels, line); // and another...

        wdCurves.emplace_back(massModel.mass);
        wdCurves.back().carbonCurves.emplace_back(newCarbon);

        while (!pCoolingModels.eof())
        {
            if (massModel.hasHLum) // Has one extra throw-away field
            {
                pCoolingModels >> ignore
                               >> newTeff
                               >> ignore >> ignore >> ignore
                               >> newAge
                               >> newRadius
                               >> ignore >> ignore >> ignore;
            }
            else
            {
                pCoolingModels >> ignore
                               >> newTeff
                               >> ignore >> ignore >> ignore
                               >> newAge
                               >> newRadius
                               >> ignore >> ignore;
            }

            if (!pCoolingModels.eof())
            {
                wdCurves.back().carbonCurves.back().records.emplace_back(newRadius, log10(1e6) + newAge, newTeff);
            }
        }
    }
}
