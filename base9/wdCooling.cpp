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

#include "evolve.hpp"
#include "linInterp.hpp"
#include "wdCooling.hpp"

using std::lower_bound;
using std::string;
using std::ifstream;
using std::vector;
using std::cerr;
using std::endl;

static vector<struct wdCoolingCurve> wdCurves;
static int coolingModel;

struct althausModel
{
    const string filename;
    const bool hasHLum;
    const double mass;
};

struct renedoModel
{
    const string filename;
    const double mass;
};

void loadWDCool (string path, int modelSet)
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

    static struct renedoModel renedo[] = {
        {"wd0524_z001.trk", 0.524},
        {"wd0570_z001.trk", 0.570},
        {"wd0593_z001.trk", 0.593},
        {"wd0609_z001.trk", 0.609},
        {"wd0632_z001.trk", 0.632},
        {"wd0659_z001.trk", 0.659},
        {"wd0705_z001.trk", 0.705},
        {"wd0767_z001.trk", 0.767},
        {"wd0837_z001.trk", 0.837},
        {"wd0877_z001.trk", 0.877},
        {"wd0934_z001.trk", 0.934},
    };

    string tempFile, line;
    double newAge, newTeff, newMass, newRadius;
    double newCarbon = 0.6 ;// 0.38; // Good default value, per Mike Montgomery
    double ignore;

    coolingModel = modelSet;

    tempFile = path;

    if (modelSet == WOOD)
    {
        tempFile += "xb.comb";
    }
    else if (modelSet == MONTGOMERY)
    {
        tempFile += "wdtables";
    }
    else if ((modelSet != ALTHAUS) && (modelSet != RENEDO))
    {
        cerr << "\nCooling models do not exist.  Exiting..." << endl;
        exit (1);
    }

    if ((modelSet == MONTGOMERY) || (modelSet == WOOD))
    {
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
            if (modelSet == WOOD)
            {
                pCoolingModels >> ignore
                               >> newAge
                               >> ignore >> ignore >> ignore
                               >> newRadius
                               >> newTeff
                               >> ignore >> ignore >> ignore
                               >> newMass;
            }
            else
            {
                pCoolingModels >> ignore
                               >> newAge
                               >> ignore >> ignore >> ignore
                               >> newRadius
                               >> newTeff
                               >> ignore >> ignore >> ignore
                               >> newMass
                               >> newCarbon;
            }

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
    else if (modelSet == ALTHAUS)
    {
        for (auto massModel : althaus)
        {
            tempFile = path;
            tempFile += "althaus/";
            tempFile += massModel.filename;

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
    else if (modelSet == RENEDO)
    {
        for (auto massModel : renedo)
        {
            tempFile = path;
            tempFile += "renedo/";
            tempFile += massModel.filename;

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
                pCoolingModels >> ignore
                               >> newTeff
                               >> ignore >> ignore
                               >> newAge
                               >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                               >> newRadius;

                if (!pCoolingModels.eof() && newAge >= 0.0)
                {
                    wdCurves.back().carbonCurves.back().records.emplace_back(log10(newRadius * R_sun), log10(1e6 + newAge), newTeff);
                }
            }
        }
    }
}


double wdMassToTeffAndRadius_montgomery (double logAge, double x_carbon, double wdPrecLogAge, double wdMass, double &thisWDLogRadius)
{
    vector<double> carbonTeff;
    vector<double> carbonRadius;

    double wdCoolLogAge = 0.0;

    if (logAge > wdPrecLogAge)
    {                           // age limit check: otherwise wdCoolLogAge
        wdCoolLogAge = log10(exp10(logAge) - exp10(wdPrecLogAge));
    }
    else
    {                           // mcmc.c can cause this by adjusting masses and ages
        thisWDLogRadius = 0.0;
        return 0.0;                     // no need to calculate anything, return to evolve.c here
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

    thisWDLogRadius = linInterpExtrap (massIter[0].mass, massIter[1].mass, carbonRadius[0], carbonRadius[1], wdMass);

    return linInterpExtrap (massIter[0].mass, massIter[1].mass, carbonTeff[0], carbonTeff[1], wdMass);
}

double wdMassToTeffAndRadius_wood (double logAge, double wdPrecLogAge, double wdMass, double &thisWDLogRadius)
{
    vector<double> ageTeff;
    vector<double> ageRadius;

    double wdCoolLogAge = 0.0;

    if (logAge > wdPrecLogAge)
    {                           // age limit check: otherwise wdCoolLogAge
        wdCoolLogAge = log10(exp10(logAge) - exp10(wdPrecLogAge));
    }
    else
    {                           // mcmc.c can cause this by adjusting masses and ages
        thisWDLogRadius = 0.0;
        return 0.0;                     // no need to calculate anything, return to evolve.c here
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
        record r(0, wdCoolLogAge, 0);
        auto ageIter = lower_bound(m->carbonCurves[0].records.begin(), m->carbonCurves[0].records.end(), r, record::compareAge);
        assert(ageIter - m->carbonCurves[0].records.begin() > 0);

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

    thisWDLogRadius = linInterpExtrap (massIter[0].mass, massIter[1].mass, ageRadius[0], ageRadius[1], wdMass);

    return linInterpExtrap (massIter[0].mass, massIter[1].mass, ageTeff[0], ageTeff[1], wdMass);
}


double wdMassToTeffAndRadius (double logAge, double x_carbon, double wdPrecLogAge, double wdMass, double &thisWDLogRadius)
{
    if ((coolingModel == WOOD) || (coolingModel == ALTHAUS) || (coolingModel == RENEDO))
    {
        return wdMassToTeffAndRadius_wood (logAge, wdPrecLogAge, wdMass, thisWDLogRadius);
    }
    else
    {
        return wdMassToTeffAndRadius_montgomery (logAge, x_carbon, wdPrecLogAge, wdMass, thisWDLogRadius);
    }
}
