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
#include "WoodWdModel.hpp"

using std::lower_bound;
using std::string;
using std::ifstream;
using std::vector;
using std::cerr;
using std::endl;

void WoodWdModel::loadModel(string path, FilterSetName)
{
    string tempFile, line;
    double newAge, newTeff, newMass, newRadius;
    double newCarbon = 0.38; // Good default value, per Mike Montgomery
    double ignore;

    tempFile = path + "xb.comb";

    ifstream pCoolingModels;
    pCoolingModels.open(tempFile);

    if (!pCoolingModels.is_open())
    {
        cerr << "\n file " << tempFile << " was not found - exiting" << endl;
        exit (1);
    }        

    getline(pCoolingModels, line); // get header line

    while (!pCoolingModels.eof())
    {
        pCoolingModels >> ignore
                       >> newAge
                       >> ignore >> ignore >> ignore
                       >> newRadius
                       >> newTeff
                       >> ignore >> ignore >> ignore
                       >> newMass;

        if (!pCoolingModels.eof())
        {
            if (wdCurves.empty() || newMass != wdCurves.back().mass)
            {
                wdCurves.emplace_back(newMass);
            }
       
            if (wdCurves.back().carbonCurves.empty() || newCarbon != wdCurves.back().carbonCurves.back().carbon)
            {
                wdCurves.back().carbonCurves.emplace_back(newCarbon);
            }

            wdCurves.back().carbonCurves.back().records.emplace_back(newRadius, log10(newAge), newTeff);
        }
    }
}
