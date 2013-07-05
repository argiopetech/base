#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cassert>

#include "evolve.hpp"
#include "linInterp.hpp"
#include "MontgomeryWdModel.hpp"

using std::lower_bound;
using std::string;
using std::ifstream;
using std::pair;
using std::make_pair;
using std::vector;
using std::cerr;
using std::endl;

void MontgomeryWdModel::loadModel(string path)
{
    string tempFile, line;
    double newAge, newTeff, newMass, newRadius;
    double newCarbon = 0.6 ;// 0.38; // Good default value, per Mike Montgomery
    double ignore;

    tempFile = path + "wdtables";

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
                       >> newMass
                       >> newCarbon;

        if (!pCoolingModels.eof())
        {
            if (wdCurves.empty() || newMass != wdCurves.back().mass)
            {
                wdCurves.emplace_back(newMass);
       
                if (wdCurves.back().carbonCurves.empty() || newCarbon != wdCurves.back().carbonCurves.back().carbon)
                {
                    wdCurves.back().carbonCurves.emplace_back(newCarbon);
                }

                wdCurves.back().carbonCurves.back().records.emplace_back(newRadius, log10(newAge), newTeff);
            }
        }
    }
}


pair<double, double> MontgomeryWdModel::wdMassToTeffAndRadius(double logAge, double x_carbon, double wdPrecLogAge, double wdMass) const
{
    vector<double> carbonTeff;
    vector<double> carbonRadius;

    double wdCoolLogAge = 0.0, thisTeff, thisRadius;

    if (logAge > wdPrecLogAge)
    {                           // age limit check: otherwise wdCoolLogAge
        wdCoolLogAge = log10(exp10(logAge) - exp10(wdPrecLogAge));
    }
    else
    {                           // mcmc.c can cause this by adjusting masses and ages
        return make_pair<double, double> (0.0, 0.0);                     // no need to calculate anything, return to evolve.c here
    }

    auto massIter = lower_bound(wdCurves.begin(), wdCurves.end(), wdCoolingCurve(wdMass));
    assert(massIter - wdCurves.begin() > 0);

    if (massIter == wdCurves.end())
    {
        massIter -= 2;
    }
    else
    {
        massIter -= 1;
    }

    //For each mass entry, interpolate in age
    for (auto m = massIter; m <= massIter + 1; ++m)
    {
        vector<double> ageTeff;
        vector<double> ageRadius;

        auto carbonIter = lower_bound(m->carbonCurves.begin(), m->carbonCurves.end(), wdCarbonCurve(x_carbon));
        assert(carbonIter - m->carbonCurves.begin() > 0);

        if (carbonIter == m->carbonCurves.end())
        {
            cerr << "Extrapolating carbon" << endl;
            carbonIter -= 2;
        }
        else
        {
            carbonIter -= 1;
        }


        for (auto c = carbonIter; c <= carbonIter + 1; ++c)
        {
            record r(0, wdCoolLogAge, 0);
            auto ageIter = lower_bound(c->records.begin(), c->records.end(), r, record::compareAge);
            assert(ageIter - c->records.begin() > 0);

            if (ageIter == m->carbonCurves[0].records.end())
            {
                cerr << "Extrapolating age" << endl;
                ageIter -= 2;
            }
            else
            {
                ageIter -= 1;
            }

            ageTeff.push_back(linInterpExtrap (ageIter[0].logAge, ageIter[1].logAge, ageIter[0].logTeff, ageIter[1].logTeff, wdCoolLogAge));

            ageRadius.push_back(linInterpExtrap (ageIter[0].logAge, ageIter[1].logAge, ageIter[0].logRadius, ageIter[1].logRadius, wdCoolLogAge));
        }

        // Now interpolate in mass
        carbonTeff.push_back(linInterpExtrap (carbonIter[0].carbon, carbonIter[1].carbon, ageTeff[0], ageTeff[1], x_carbon));

        carbonRadius.push_back(linInterpExtrap (carbonIter[0].carbon, carbonIter[1].carbon, ageRadius[0], ageRadius[1], x_carbon));
    }

    thisRadius = linInterpExtrap (massIter[0].mass, massIter[1].mass, carbonRadius[0], carbonRadius[1], wdMass);
    thisTeff = linInterpExtrap (massIter[0].mass, massIter[1].mass, carbonTeff[0], carbonTeff[1], wdMass);

    return pair<double, double> (thisTeff, thisRadius);
}
