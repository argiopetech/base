#include <algorithm>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <cmath>

#include "Cluster.hpp"
#include "LinearTransform.hpp"
#include "Matrix.hpp"
#include "Star.hpp"
#include "YaleMsModel.hpp"

using std::array;
using std::getline;
using std::ifstream;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::vector;

const unsigned int maxIgnore = std::numeric_limits<char>::max();

//global variables
// static int iZ, iAge;
// static double yyFeH[N_YY_Z], yyZ[N_YY_Z];
// static double yyLogAge[N_YY_Z][N_YY_AGES], yyAge[N_YY_Z][N_YY_AGES];
// static struct yyIsochrone yyIso[N_YY_Z][N_YY_AGES];
// static double yyAGBt[N_YY_Z][N_YY_AGES];
// static struct globalIso tempIso[2] {globalIso(MAX_YY_ENTRIES, N_YY_FILTS), globalIso(MAX_YY_ENTRIES, N_YY_FILTS)};
static Matrix<double, 6, 8> coeff;

//Static funtions
static void eepset (double x[], double y[], int nstep, int igrd, double *slope);
static int ToffM (double x[4], double yy[4], double *xp, double *ybar, int iorder);
static double feh2z (double FeH);


static vector<string> getFileNames (const string path)
{
    const array<string, 11> fileNames = {
        "76997z00001a0o2v2",
        "7697z0001a0o2v2",
        "7688z0004a0o2v2",
        "767z001a0o2v2",
        "758z004a0o2v2",
        "749z007a0o2v2",
        "74z01a0o2v2",
        "71z02a0o2v2",
        "65z04a0o2v2",
        "59z06a0o2v2",
        "53z08a0o2v2"
    };

    vector<string> files;

    for (auto file : fileNames)
    {
        string tempFile = path;

        tempFile += "YYiso/yy00l.x";

        tempFile += file;

        files.push_back(tempFile);
    }

    return files;

}


string YaleMsModel::getFileName (string) const
{
    throw std::logic_error("YY models do not support single-file loading");
}


bool YaleMsModel::isSupported(FilterSetName filterSet) const
{
    return (filterSet == FilterSetName::UBVRIJHK);
}


void YaleMsModel::loadModel (const string path, FilterSetName filterSet)
{
    assert(isSupported(filterSet));

    ifstream fin;

    for (auto file : getFileNames(path))
    {
        double fileZ, fileFeH, logAge, ignore;

        string line;
        stringstream lin;

        vector<Isochrone> isochrones;
        vector<EvolutionaryPoint> eeps;

        // Open the file
        fin.open(file);

        // Ensure the file opened successfully
        if (!fin)
        {
            cerr << "\n file " << file << " was not found - exiting" << endl;
            exit (1);
        }

        // Get Z from the first header line
        fin.ignore(maxIgnore, '='); // Eat the first bit of the line
        fin >> fileZ;

        // Get FeH from the first header line
        // Z=0.080000 Y=0.390000 OS=0.20 l/Hp=1.743201 [Fe/H]= 0.305363 [Alpha/Fe]= 0.60
        fin.ignore(maxIgnore, '='); // Y value
        fin.ignore(maxIgnore, '='); // OS value
        fin.ignore(maxIgnore, '='); // l/Hp value
        fin.ignore(maxIgnore, '='); // [Fe/H] value
        fin >> fileFeH;

        fin.ignore(maxIgnore, '\n'); // Eat the last bit of the line
        getline(fin, line); // Eat the second header line (column headers)

        // Now, for every other line in the file...
        while (!fin.eof())
        {
            getline(fin, line);

            if ((line.size() > 1) && (line.at(0) != 'a'))
            {
                double tempMass;

                stringstream in(line);

                array<double, FILTS> mags;
                array<double, 7> params;
                mags.fill(99.999);

                in >> tempMass
                   >> ignore >> ignore >> ignore
                   >> mags.at(2); // V = Mv

                // U-B    B-V    V-R    V-I    V-J    V-H    V-K
                for (auto &param : params)
                {
                    in >> param;
                }

                mags.at(1) = params.at(1) + mags.at(2); // B = (B-V) + V
                mags.at(0) = params.at(0) + mags.at(1); // U = (U-B) + B

                for (int i = 2; i < 7; ++i)
                {
                    mags.at(i + 1) =  mags.at(2) - params.at(i); // R = V - (V-R), etc
                }

                if (!fin.eof())
                {
                    // Yale doesn't have EEPs, so we'll pretend with 0
                    eeps.emplace_back(0, tempMass, mags);
                }
            }
            else if (line.size() == 1) // End of an isochrone
            {
                isochrones.emplace_back(logAge, eeps);
                eeps.clear();
            }
            else if (line.at(0) == 'a') // Beginning of a new isochrone
            {
                stringstream in(line);

                in.ignore(maxIgnore, '='); // Eat the first bit of the line
                in >> logAge;

                logAge = log10(logAge * 1e9); // logAge is actually in non-log gigayears in the file
            }
        } // EOF

//        fileZ = log10(fileZ / modelZSolar);
        vector<HeliumCurve> tVector;
        tVector.emplace_back(0.0, isochrones);
        fehCurves.emplace_back(fileFeH, tVector); // Push the entire isochrone set into the model's FehCurve vector

        fin.close();
    }

    /////////////////////////////////////////////////////////////////////
    // Read in the coefficients needed to calculate the precurser ages //
    /////////////////////////////////////////////////////////////////////

    // Open coeff file for reading
    auto tempFile = path + "YYiso/yyAGBtcoeff.dat";
    fin.open(tempFile);

    //fscanf(pModelList,"%s",tempFile);
    if (!fin)
    {
        cerr << "\n file " << tempFile << " was not found - exiting" << endl;
        exit (1);
    }

    // Read in the coefficients needed to calculate the precurser ages
    for (auto &i : coeff)
    {
        for (auto &j : i)
        {
            fin >> j;
        }
    }

    fin.close();

    //Set the min and max age in this model set (for use in densities.c)
    ageLimit.first = fehCurves.front().heliumCurves.front().isochrones.front().logAge;
    ageLimit.second = fehCurves.front().heliumCurves.front().isochrones.back().logAge;

    assert(ageLimit.second == fehCurves.back().heliumCurves.back().isochrones.back().logAge);
}


double YaleMsModel::deriveAgbTipMass (const vector<int> &filters, double newFeH, double ignored, double newAge)
{
    isochrone = deriveIsochrone(filters, newFeH, ignored, newAge);

    return isochrone.agbTipMass();
}


Isochrone YaleMsModel::deriveIsochrone (const vector<int> &filters, double newFeH, double, double newAge) const
{
    if ((newAge < fehCurves.front().heliumCurves.front().isochrones.front().logAge)
     || (newAge > fehCurves.back().heliumCurves.front().isochrones.back().logAge)
     || (newFeH < fehCurves.front().feh)
     || (newFeH > fehCurves.back().feh))
    {
        throw InvalidCluster("Age or FeH out of bounds in YaleMsModel::deriveIsochrone");
    }


//     double newAge = exp10 (newLogAge) / 1e9;
//     double newZ;

//     newZ = feh2z (newFeH);

//     iAge = -1;
//     iZ = -1;

//     // Find the values for each parameter that we will be interpolating
//     // between and calculate the interpolation coefficients.
//     iZ = binarySearch (yyZ, N_YY_Z, newZ);
//     iAge = binarySearch (yyAge[iZ], N_YY_AGES, newAge);

//     iZ--;
//     if (iZ >= N_YY_Z - 4)
//         iZ = N_YY_Z - 4;
//     if (iZ <= 0)
//         iZ = 0;
// //    newY = dydz*(newZ-zp)+yp;

//     //Will eventually interpolate in age between these two isochrones
//     for (int i = 0; i < 2; i++)
//     {
//         for (int m = 0; m < yyIso[iZ + 1][i + iAge].nEntries; m++)
//         {
//             tempIso[i].mass[m] = CUBEINT (log (yyZ[iZ]), yyIso[iZ][i + iAge].mass[m], log (yyZ[iZ + 1]), yyIso[iZ + 1][i + iAge].mass[m], log (yyZ[iZ + 2]), yyIso[iZ + 2][i + iAge].mass[m], log (yyZ[iZ + 3]), yyIso[iZ + 3][i + iAge].mass[m], log (newZ));

//             for (int p = 0; p < N_YY_FILTS; p++)
//             {
//                 tempIso[i].mag[m][p] = CUBEINT (log (yyZ[iZ]), yyIso[iZ][i + iAge].mag[m][p], log (yyZ[iZ + 1]), yyIso[iZ + 1][i + iAge].mag[m][p], log (yyZ[iZ + 2]), yyIso[iZ + 2][i + iAge].mag[m][p], log (yyZ[iZ + 3]), yyIso[iZ + 3][i + iAge].mag[m][p], log (newZ));
//             }
//         }

//         tempIso[i].nEntries = yyIso[iZ + 1][iAge + i].nEntries;
//         tempIso[i].age = yyIso[iZ][iAge + i].age;
//         tempIso[i].logAge = yyIso[iZ][iAge + i].logAge;
//     }

//     //SD -- Isochrones for the lower ages have fewer entries (~35 instead of 140)
//     //If this age is right on the border, just use the two higher age entries
//     //that have all 140 entries (and extrapolate, technically)
//     int iYYm = tempIso[0].nEntries; //yyIso[iZ][iAge+1].nEntries;
//     if (yyIso[iZ][iAge].nEntries < iYYm)
//         iAge++;

//     for (int m = 0; m < iYYm; m++)
//     {
//         isochrone.mass[m] = linearTransform<>(tempIso[0].age, tempIso[1].age, tempIso[0].mass[m], tempIso[1].mass[m], newAge).val;

//         for (int p = 0; p < N_YY_FILTS; p++)
//         {
//             isochrone.mag[m][p] = linearTransform<>(tempIso[0].age, tempIso[1].age, tempIso[0].mag[m][p], tempIso[1].mag[m][p], newAge);
//         }
//     }

//     isochrone.nEntries = iYYm;

//     isochrone.age = newAge;
//     isochrone.logAge = log10 (newAge * 1e9);
//     isochrone.z = newZ;
//     isochrone.AgbTurnoffMass = isochrone.mass[isochrone.nEntries - 1];

//     assert(is_sorted(isochrone.mass.begin(), isochrone.mass.begin() + isochrone.nEntries));

    return {};
}


static void eepset (double x[], double y[], int nstep, int igrd, double *slope)
{
    // double xp = 0.0, ybar = 0.0, ymax = 0.0;
    // double tomass, /*totemp, */ anchor, xhi, xlo, eeptag;
    // int ikip = 1, i = 0;
    // int rv = -1;

    // //To find the turnoff mass
    // for (i = 0; i < nstep; i++)
    // {
    //     if (y[i] >= ymax)
    //     {
    //         ikip = i;
    //         ymax = y[i];
    //     }
    // }

    // if (ToffM (&(x[ikip - 1]), &(y[ikip - 1]), &xp, &ybar, 0))
    // {
    //     xp = x[ikip];
    // }
    // else
    // {
    //     rv = 0;
    //     if (xp < x[ikip])
    //     {
    //         rv = ToffM (&(x[ikip - 2]), &(y[ikip - 2]), &xp, &ybar, 1);
    //     }
    //     if (rv)
    //     {
    //         xp = x[ikip];
    //     }
    //     else
    //     {
    //         if (xp < x[ikip] && xp > x[ikip + 1])
    //             cerr << "possibly incorrect" << endl;
    //     }
    // }

    // tomass = xp;

    // // Got the turnoff mass
    // // To find the slope
    // anchor = (tomass - x[0]) / (x[nstep - 1] - x[0]);
    // xhi = 5.0;
    // xlo = 0.005;
    // (*slope) = 0.5 * (xhi + xlo);

    // while (1)
    // {
    //     eeptag = atan ((*slope) * (igrd - 1)) / atan ((*slope) * (nstep - 1));
    //     if (fabs (anchor - eeptag) <= 0.5e-7)
    //         break;
    //     if (eeptag >= anchor)
    //         xhi = (*slope);
    //     else
    //         xlo = (*slope);
    //     (*slope) = 0.5 * (xhi + xlo);
    //     if (fabs ((xhi - xlo) / (*slope)) <= 0.5e-12)
    //         break;
    // }

    // // This is never read, so we're getting rid of it. Probably not a good thing.
    // // eeptag = atan((*slope)*(igrd-1))/atan((*slope)*(nstep-1));

    // return;
}

//SD not to self, makesure iorder is fed in as one less than in Fortran code
static int ToffM (double x[4], double yy[4], double *xp, double *ybar, int iorder)
{
    // int n = 4, k, l, m;
    // double a, b, c, s, xm, xbar;

    // vector<double> y(n);

    // for (k = 0; k < n; k++)
    //     y.at(k) = yy[k];

    // for (k = 0; k < n - 1; k++)
    // {
    //     for (l = 0; l < n - k - 1; l++)
    //     {
    //         y.at(l) = (y.at(l + 1) - y.at(l)) / (x[l + k + 1] - x[l]);
    //     }
    // }

    // a = 3.0 * y.at(0);
    // b = -2.0 * y.at(0) * (x[3] + x[2] + x[1]) + 2.0 * y.at(1);
    // c = y.at(0) * (x[2] * x[3] + x[1] * x[3] + x[1] * x[2]) - y.at(1) * (x[3] + x[2]) + y.at(2);
    // s = sqrt (b * b - 4.0 * a * c);

    // (*xp) = (-b + s) / (2.0 * a);
    // xm = (-b - s) / (2.0 * a);

    // if ((*xp) >= x[iorder] && (*xp) <= x[2])
    //     xbar = (*xp);
    // else if (xm >= x[iorder] && xm <= x[2])
    //     xbar = xm;
    // else
    // {
    //     if (iorder >= 1)
    //         return 1;
    //     else
    //     {
    //         cerr << "Failed to find the turnoff mass (xp=" << *xp << ", xm=" << xm << ").  Exiting." << endl;
    //         exit (1);
    //     }
    // }

    // (*ybar) = y.at(0);

    // for (m = 1; m < n; m++)
    // {
    //     (*ybar) = (*ybar) * (xbar - x[m]) + y.at(m);
    // }

    // (*xp) = xbar;

    return 0;
}

static double feh2z (double FeH)
{
    double Z = 0.0, FeH0 = 0.0, ZovX = 0.0;
    double YYaf = 0.0;

    FeH0 = FeH - QUAD (0., 0., 0.3, FeHa2, 0.6, FeHa4, YYaf);
    ZovX = exp10(FeH0) * Zsun / Xsun;
    Z = ZovX * (1.0 + dydz * zp - yp) / (1.0 + ZovX * (1.0 + dydz));

    return Z;
}


/*************************************************************************************
Determine WD precursor age by using coefficients derived from fitting to the current
Yale-Yonsei models.  If the models change, these will need to change as well.
Note that the appropriate AgbTurnoffMass mass and lifetime is not the ZAMS mass and lifetime of
the star currently at the AgbTurnoffMass, but rather refers to the properties of the potentially
higher mass and younger AgbTurnoffMass star that was the WD precursor.
*************************************************************************************/
double YaleMsModel::wdPrecLogAge (double thisFeH, double zamsMass)
{
    if (zamsMass < 1.0)
    {
        zamsMass = 1.0;
    }
    else if (zamsMass > 8.0)
    {
        zamsMass = 8.0;
    }

    double wdPrecLogAge = 0.0;

    for (int i = 0; i < 6; i++)
    {
        double maCoeff = 0.0;

        for (int j = 0; j < 8; j++)
        {
            maCoeff += coeff[i][j] * pow (thisFeH, j);
        }

        wdPrecLogAge += maCoeff * pow (zamsMass, i);
    }

    return wdPrecLogAge;
}
