#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <unistd.h>

#include "evolve.hpp"
#include "loadModels.hpp"
#include "Settings.hpp"
#include "FilterSet.hpp"

using std::vector;

const int CLUS_READ       =  9;
const int CLUS_STAT_WRITE = 10;
const int MASS1_FILE      = 11;
const int MASS2_FILE      = 12;
const int CMD_FILE        = 13;
const int DEBUG_FILE      = 16;

// Used by evolve.c
double ltau[2], wdLogTeff[2];
int aFilt = 0;

extern struct globalIso isochrone;

// Used by a bunch of different functions.
int verbose, needMassNow = 1, useFilt[FILTS], numFilts;

static void openOutputFiles (FILE ** filePtr, char *filename, int fileType);

int main (int argc, char *argv[])
{
    int j, filt, iMag;
    double minMag, maxMag;
    char line[240], filename[100];

    FILE *rDataPtr;
    FILE *wDebugPtr;

    Cluster theCluster;
    vector<Star> stars;

    Settings settings;

    settings.fromCLI (argc, argv);
    if (!settings.files.config.empty())
    {
        settings.fromYaml (settings.files.config);
    }
    else
    {
        settings.fromYaml ("base9.yaml");
    }

    settings.fromCLI (argc, argv);

    Model evoModels = makeModel(settings);

    ////////////////////////////////////////////
    /////// Open files to read and write ///////
    ////////////////////////////////////////////

    /* printf("Enter file name containing color/magnitude data:\n> "); */
    /* scanf("%s",filename); */
    strcpy (filename, settings.files.phot.c_str());
    if ((rDataPtr = fopen (filename, "r")) == NULL)
    {
        printf ("***Error: file %s was not found.***\n", filename);
        printf ("[Exiting...]\n");
        exit (1);
    }

    /* printf("Enter minimum and maximum magnitude of MS to use and band (0,1, or 2):\n> "); */
    /* scanf("%lf %lf %d", &minMag, &maxMag,&iMag); */
    minMag = settings.cluster.minMag;
    maxMag = settings.cluster.maxMag;
    iMag = settings.cluster.index;
    if (iMag < 0 || iMag > 2)
    {
        printf ("***Error: %d not a valid magnitude index.  Choose 0,1,or 2.***\n", iMag);
        printf ("[Exiting...]\n");
        exit (1);
    }

    // char outfilename[100];
    /* printf("\n Enter isochrone file name : "); */
    /* scanf("%s",filename); */
    strcpy (filename, settings.files.output.c_str());  // File ending gets added in openOutputFiles
    openOutputFiles (&wDebugPtr, filename, DEBUG_FILE);

    ////////////////////////////////////////////
    ////// Determine which filters to use //////
    ////// read header from .scatter file //////
    ////////////////////////////////////////////

    fgets (line, 240, rDataPtr);        /* skip first header line */
    for (filt = 0; filt < FILTS; filt++)
    {
        fscanf (rDataPtr, "%d ", &useFilt[filt]);
        if (useFilt[filt])
        {
            evoModels.numFilts++;
            //printf("** %d **\n",evoModels.numFilts);
            aFilt = filt;               // Sets this to a band we know we are using (for evolve)
        }
    }
    fgets (line, 240, rDataPtr);        // and next header line
    fclose (rDataPtr);

    ///////////////////////////////////
    ///// Read in other variables /////
    ///////// and load models /////////
    ///////////////////////////////////
    theCluster.M_wd_up = settings.whiteDwarf.M_wd_up;

    verbose = settings.verbose;
    if (verbose < 0 || verbose > 2)
        verbose = 1;            /* give standard feedback if incorrectly specified */


    loadModels (&theCluster, evoModels, settings);      /* read in stellar evol & WD models */


    /* read cluster parameters */
    theCluster.parameter[AGE] = settings.cluster.logClusAge;

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

    // Create a simulated cluster based on the cluster stats and output for comparison
    theCluster.nStars = 1;

    stars.resize(theCluster.nStars);

    stars.at(0).U = 1.0;
    stars.at(0).massRatio = 0.0;

    evolve (&theCluster, evoModels, stars, -1);

    int nWD = 10000;

    theCluster.nStars = isochrone.nEntries + nWD;

    stars.resize(theCluster.nStars);

    double dMass = (theCluster.M_wd_up - isochrone.mass[isochrone.nEntries]) / (double) nWD;

    for (j = 0; j < theCluster.nStars; j++)
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
    evolve (&theCluster, evoModels, stars, -1);


    fprintf (wDebugPtr, " mass stage1");
    for (filt = 0; filt < FILTS; filt++)
        if (useFilt[filt])
            fprintf (wDebugPtr, "          %s", getFilterName (filt));
    fprintf (wDebugPtr, "\n");

//   double prevV=100;
    for (j = 0; j < theCluster.nStars; j++)
    {

        if (stars.at(j).photometry[0] < 90)
        {
            if (stars.at(j).status[0] == MSRG)
            {
                // if(prevV > stars.at(j).photometry[0]){
                fprintf (wDebugPtr, "%lf %6d ", stars.at(j).U, stars.at(j).status[0]);
                for (filt = 0; filt < evoModels.numFilts; filt++)
                    fprintf (wDebugPtr, "%10f ", stars.at(j).photometry[filt]);
                fprintf (wDebugPtr, "\n");
//           prevV = stars.at(j).photometry[0];
                // }
                // else{
                //   prevV = 0.0;
                //   continue;
                // }
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
            printf ("***Error: Bad file choice in openOutputFiles().***\n");
            printf ("[Exiting...]\n");
            exit (0);

    }

    printf ("Reading file  : %s (%s)\n", tmpfile, mode);
    if ((*filePtr = fopen (tmpfile, mode)) == NULL)
    {
        printf ("***Error: File %s was not available. %s ***\n", tmpfile, mode);
        printf ("[Exiting...]\n");
        exit (1);
    }
}
