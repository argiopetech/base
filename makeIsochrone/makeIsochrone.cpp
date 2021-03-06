#include <array>
#include <iostream>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <unistd.h>

#include "evolve.hpp"
#include "Settings.hpp"
#include "FilterSet.hpp"

using std::array;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

const int CLUS_READ       =  9;
const int CLUS_STAT_WRITE = 10;
const int MASS1_FILE      = 11;
const int MASS2_FILE      = 12;
const int CMD_FILE        = 13;
const int DEBUG_FILE      = 16;

vector<int> filters;

// Used by evolve.c
double wdLogTeff[2];

double filterPriorMin[FILTS];
double filterPriorMax[FILTS];


// Used by a bunch of different functions.
static void openOutputFiles (FILE ** filePtr, char *filename, int fileType);

int main (int argc, char *argv[])
{
    int filt, iMag;
    double minMag, maxMag;
    char line[240], filename[100];

    FILE *rDataPtr;
    FILE *wDebugPtr;

    array<double, 2> ltau;

    Cluster theCluster;
    vector<Star> stars;
    array<double, FILTS> globalMags;

    struct globalIso isochrone;

    Settings settings;

    settings.loadSettings (argc, argv);

    if (settings.seed == std::numeric_limits<uint32_t>::max())
    {
        srand(std::time(0));
        settings.seed = rand();

        cout << "Seed: " << settings.seed << endl;
    }

    const Model evoModels = makeModel(settings);

    ////////////////////////////////////////////
    /////// Open files to read and write ///////
    ////////////////////////////////////////////

    strcpy (filename, settings.files.phot.c_str());
    if ((rDataPtr = fopen (filename, "r")) == NULL)
    {
        cerr << "***Error: file " << filename << " was not found.***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }

    minMag = settings.cluster.minMag;
    maxMag = settings.cluster.maxMag;
    iMag = settings.cluster.index;
    if (iMag < 0 || iMag > 2)
    {
        cerr << "***Error: " << iMag << " not a valid magnitude index.  Choose 0,1,or 2.***" << endl;
    cerr << "[Exiting...]" << endl;
        exit (1);
    }

    strcpy (filename, settings.files.output.c_str());  // File ending gets added in openOutputFiles
    openOutputFiles (&wDebugPtr, filename, DEBUG_FILE);

    ////////////////////////////////////////////
    ////// Determine which filters to use //////
    ////// read header from .scatter file //////
    ////////////////////////////////////////////

    fgets (line, 240, rDataPtr);        /* skip first header line */
    for (filt = 0; filt < FILTS; filt++)
    {
        int useFilt = false;
        fscanf (rDataPtr, "%d ", &useFilt);
        if (useFilt)
        {
            filters.push_back(filt);
            const_cast<Model&>(evoModels).numFilts++;
        }
    }
    fgets (line, 240, rDataPtr);        // and next header line
    fclose (rDataPtr);

    ///////////////////////////////////
    ///// Read in other variables /////
    ///////// and load models /////////
    ///////////////////////////////////
    theCluster.M_wd_up = settings.whiteDwarf.M_wd_up;
    theCluster.carbonicity = settings.whiteDwarf.carbonicity;

    /* read cluster parameters */
    theCluster.parameter[AGE] = settings.cluster.logAge;

    theCluster.parameter[FEH] = settings.cluster.Fe_H;

    theCluster.parameter[MOD] = settings.cluster.distMod;

    theCluster.parameter[ABS] = settings.cluster.Av;

    if (settings.mainSequence.msRgbModel == MsModel::CHABHELIUM)
    {
        theCluster.parameter[YYY] = settings.cluster.Y;
    }
    else
    {
        theCluster.parameter[YYY] = 0.0;
    }

    stars.resize(1);

    stars.front().U = 1.0;
    stars.front().massRatio = 0.0;

    for (auto s : stars)
        evolve (theCluster, evoModels, globalMags, filters, s, ltau);

    int nWD = 10000;

    stars.resize(isochrone.nEntries + nWD);

    double dMass = (theCluster.M_wd_up - isochrone.mass[isochrone.nEntries]) / (double) nWD;

    for (decltype(stars.size()) j = 0; j < stars.size(); j++)
    {
        if (j < isochrone.nEntries)
        {
            stars.at(j).U = isochrone.mass[j];
        }
        else
        {
            stars.at(j).U = stars.at(j - 1).U + dMass;
        }
        stars.at(j).massRatio = 0.0;
    }

    for (auto s : stars)
        evolve (theCluster, evoModels, globalMags, filters, s, ltau);


    fprintf (wDebugPtr, " mass stage1");
    for (auto f : filters)
        fprintf (wDebugPtr, "          %s", evoModels.filterSet->getFilterName(f).c_str());
    fprintf (wDebugPtr, "\n");

    for (decltype(stars.size()) j = 0; j < stars.size(); j++)
    {

        if (stars.at(j).photometry[0] < 90)
        {
            if (stars.at(j).status[0] == MSRG)
            {
                fprintf (wDebugPtr, "%lf %6d ", stars.at(j).U, stars.at(j).status[0]);
                for (filt = 0; filt < evoModels.numFilts; filt++)
                    fprintf (wDebugPtr, "%10f ", stars.at(j).photometry[filt]);
                fprintf (wDebugPtr, "\n");
            }
            else
            {
                fprintf (wDebugPtr, "%5.2f %3d ", stars.at(j).U, stars.at(j).status[0]);
                for (filt = 0; filt < evoModels.numFilts; filt++)
                    fprintf (wDebugPtr, "%10f ", stars.at(j).photometry[filt]);
                fprintf (wDebugPtr, "\n");
            }
        }
    }

    fclose (wDebugPtr);

    return (0);
}


static void openOutputFiles (FILE ** filePtr, char *filename, int fileType)
{
    char tmpfile[100];
    const char *mode = "r";

    strcpy (tmpfile, filename);
    switch (fileType)
    {
        // output files differentiated by file name extension
        case CLUS_READ:
            strcat (tmpfile, ".cluster");
            break;
        case MASS1_FILE:
            strcat (tmpfile, ".mass1");
            break;
        case MASS2_FILE:
            strcat (tmpfile, ".mass2");
            break;
        case CLUS_STAT_WRITE:
            strcat (tmpfile, ".cluster.stat");
            mode = "w";
            break;
        case CMD_FILE:
            strcat (tmpfile, ".cmd");
            mode = "w";
            break;
        case DEBUG_FILE:
            strcat (tmpfile, ".cmd.debug");
            mode = "w";
            break;
        default:
            cerr << "***Error: Bad file choice in openOutputFiles().***" << endl;
            cerr << "[Exiting...]" << endl;
            exit (0);

    }

    cout << "Reading file  : " << tmpfile << " (" << mode << ")" << endl;
    if ((*filePtr = fopen (tmpfile, mode)) == NULL)
    {
        cerr << "***Error: File " << tmpfile << " was not available. " << mode << "***" << endl;
        cerr << "[Exiting...]" << endl;
        exit (1);
    }
}
