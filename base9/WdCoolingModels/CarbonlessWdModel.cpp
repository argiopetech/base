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
#include "CarbonlessWdModel.hpp"

using std::lower_bound;
using std::string;
using std::ifstream;
using std::pair;
using std::make_pair;
using std::vector;
using std::cerr;
using std::endl;

pair<double, double> CarbonlessWdModel::wdMassToTeffAndRadius (double logAge, double, double wdPrecLogAge, double wdMass) const
{
    vector<double> ageTeff;
    vector<double> ageRadius;

    double wdCoolLogAge = 0.0, thisTeff, thisRadius;

    if (logAge > wdPrecLogAge)
    {                           // age limit check: otherwise wdCoolLogAge
        wdCoolLogAge = log10(exp10(logAge) - exp10(wdPrecLogAge));
    }
    else
    {                           // mcmc.c can cause this by adjusting masses and ages
        return make_pair<double, double> (0.0, 0.0); // no need to calculate anything, return to evolve.c here
    }

    auto massIter = lower_bound(wdCurves.begin(), wdCurves.end(), wdCoolingCurve(wdMass));
    
    if (massIter == wdCurves.begin())
    {
        // log << "Mass underflow in WD Model" << endl;
    }
    else if (massIter == wdCurves.end())
    {
        // log << "Mass overflow in WD Model" << endl;
        massIter -= 2;
    }
    else
    {
        massIter -= 1;
    }

    //For each mass entry, interpolate in age
    for (auto m = massIter; m <= massIter + 1; ++m)
    {
        record r(0, wdCoolLogAge, 0);
        auto ageIter = lower_bound(m->carbonCurves[0].records.begin(), m->carbonCurves[0].records.end(), r, record::compareAge);


        if (ageIter == m->carbonCurves[0].records.begin())
        {
            // log << "Age underflow in WD Model" << endl;
        }
        else if (ageIter == m->carbonCurves[0].records.end())
        {
            // log << "Age overflow in WD Model" << endl;
            ageIter -= 2;
        }
        else
        {
            ageIter -= 1;
        }

        ageTeff.push_back(linInterpExtrap (ageIter[0].logAge, ageIter[1].logAge, ageIter[0].logTeff, ageIter[1].logTeff, wdCoolLogAge));

        ageRadius.push_back(linInterpExtrap (ageIter[0].logAge, ageIter[1].logAge, ageIter[0].logRadius, ageIter[1].logRadius, wdCoolLogAge));
    }

    thisRadius = linInterpExtrap (massIter[0].mass, massIter[1].mass, ageRadius[0], ageRadius[1], wdMass);
    thisTeff = linInterpExtrap (massIter[0].mass, massIter[1].mass, ageTeff[0], ageTeff[1], wdMass);

    return pair<double, double>(thisTeff, thisRadius);
}
