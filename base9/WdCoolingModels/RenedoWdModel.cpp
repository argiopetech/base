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

#include "RenedoWdModel.hpp"

using std::lower_bound;
using std::string;
using std::ifstream;
using std::vector;
using std::cerr;
using std::endl;

struct renedoModel
{
    const string filename;
    const double mass;
};

void RenedoWdModel::loadModel (string path)
{
    static struct renedoModel renedo[] = {
        {"wd0524_z001.trk", 0.52490},
        {"wd0570_z001.trk", 0.57015},
        {"wd0593_z001.trk", 0.59316},
        {"wd0609_z001.trk", 0.60959},
        {"wd0632_z001.trk", 0.63229},
        {"wd0659_z001.trk", 0.65988},
        {"wd0705_z001.trk", 0.70511},
        {"wd0767_z001.trk", 0.76703},
        {"wd0837_z001.trk", 0.83731},
        {"wd0877_z001.trk", 0.87790},
        {"wd0934_z001.trk", 0.93477},
    };

    string tempFile, line;
    double newAge, newTeff, newRadius;
    double newCarbon = 0.38; // Good default value, per Mike Montgomery
    double ignore;

    for (auto massModel : renedo)
    {
        tempFile = path + "renedo/" + massModel.filename;

        ifstream pCoolingModels;
        pCoolingModels.open(tempFile);

        if (!pCoolingModels.is_open())
        {
            cerr << "\n file " << tempFile << " was not found - exiting" << endl;
            exit (1);
        }        

        getline(pCoolingModels, line); // get header line

        wdCurves.emplace_back(massModel.mass);
        wdCurves.back().carbonCurves.emplace_back(newCarbon);

        while (!pCoolingModels.eof())
        {
            pCoolingModels >> ignore
                           >> newTeff
                           >> ignore >> ignore
                           >> newAge
                           >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                           >> newRadius;

            if (!pCoolingModels.eof() && newAge >= 0.0)
            {
                wdCurves.back().carbonCurves.back().records.emplace_back(log10(newRadius * R_sun), log10(1e6 * newAge), newTeff);
            }
        }
    }
}
